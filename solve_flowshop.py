"""Permutation flow-shop makespan MILP (Bowman / DP-based model) solved with Gurobi.

The model places every job on every machine in discrete time slots using binary
time-index variables ``y[i, k, t]`` (job ``i`` occupies machine ``k`` at time
``t``), plus precedence and sequencing constraints to enforce a permutation
flow-shop schedule. The objective is to minimise the makespan ``C_max``.

This is a single reusable implementation of the same model that was coded
(cell by cell) in the original assignment notebook, looped over a set of
VFR benchmark instances rather than duplicated once per file.

Coursework: Optimization Implementation in Production and Logistics (OVGU
Magdeburg), Assignment 2 — flow-shop MILP with Gurobi.
"""

import glob
import os

import gurobipy as gp
from gurobipy import GRB
import numpy as np


def read_instance(filename: str) -> tuple[int, int, np.ndarray]:
    """Read a VFR flow-shop instance file.

    Format:
        Line 1:  ``[n_jobs] [n_machines]``
        Then, per job, n_machines ``machine time`` pairs (machines 0-based).

    Returns:
        ``(n_jobs, n_machines, processing_times)`` where
        ``processing_times[i, k]`` is the processing time of job ``i`` on
        machine ``k``.
    """
    with open(filename) as f:
        lines = f.readlines()
        n, m = map(int, lines[0].strip().split())
        processing_times = np.zeros((n, m), dtype=int)
        for i in range(n):
            data = list(map(int, lines[i + 1].strip().split()))
            for j in range(0, len(data), 2):
                machine = data[j]
                time = data[j + 1]
                processing_times[i][machine] = time
    return n, m, processing_times


def build_flowshop_model(
    n: int, m: int, p_ik: np.ndarray, time_limit: float = 900.0, mip_gap: float = 0.20
) -> gp.Model:
    """Build the permutation flow-shop makespan MILP.

    Decision variables:
        ``y[i, k, t]``        binary   job i occupies machine k at time t
        ``C_ik[i, k]``        continuous completion time of job i on machine k
        ``z_ij[i, j]``        binary   job i is processed before job j on all
                                       machines (sequencing)
        ``C_max``             continuous makespan

    Constraints:
        1. machine capacity  — at most one job per machine per time slot
        2. processing demand — job i must run p_ik[i,k] slots on machine k
        3. precedence        — job i cannot start on machine k+1 before it has
                               run on machine k
        4. contiguity        — the processing blocks of a job are contiguous
        5. completion        — the last occupied slot of job i on machine k
                               defines its completion time
        6.-7. sequencing     — for every pair of jobs, one precedes the other
                               on each machine (Big-M enforced)
        8. makespan          — C_max >= completion time of every last-machine job

    Returns:
        An un-optimised Gurobi model.
    """
    # Time horizon / Big-M: any schedule finishes within n * max(p_ik).
    T = int(n * np.max(p_ik))
    M = T

    model = gp.Model("Bowman_Permutation_Flow_Shop")

    y = model.addVars(n, m, T, vtype=GRB.BINARY, name="y")
    C_ik = model.addVars(n, m, vtype=GRB.CONTINUOUS, name="C_ik")
    z_ij = model.addVars(n, n, vtype=GRB.BINARY, name="z")
    C_max = model.addVar(vtype=GRB.CONTINUOUS, name="C_max")

    model.setObjective(C_max, GRB.MINIMIZE)

    # 1. Machine capacity: only one job at a time on each machine.
    for k in range(m):
        for t in range(T):
            model.addConstr(
                gp.quicksum(y[i, k, t] for i in range(n)) <= 1,
                name=f"Capacity_m{k}_t{t}",
            )

    # 2. Processing: job i must run at least p_ik[i,k] slots on machine k.
    for i in range(n):
        for k in range(m):
            model.addConstr(
                gp.quicksum(y[i, k, t] for t in range(T)) >= p_ik[i, k],
                name=f"Processing_j{i}_m{k}",
            )

    # 3. Precedence: machine k+1 only after machine k, per job.
    for i in range(n):
        for k in range(m - 1):
            for t in range(T):
                model.addConstr(
                    p_ik[i, k] * y[i, k + 1, t]
                    <= gp.quicksum(y[i, k, l] for l in range(t)),
                    name=f"Precedence_j{i}_m{k}_t{t}",
                )

    # 4. Contiguous processing blocks.
    for i in range(n):
        for k in range(m):
            for t in range(T - 1):
                sum_after = gp.quicksum(y[i, k, l] for l in range(t + 2, T))
                model.addConstr(
                    p_ik[i, k] * (y[i, k, t] - y[i, k, t + 1]) + sum_after
                    <= p_ik[i, k],
                    name=f"Contiguous_j{i}_m{k}_t{t}",
                )

    # 5. Completion time: the occupancy at time t implies C_ik >= t+1.
    for i in range(n):
        for k in range(m):
            for t in range(T):
                model.addConstr(
                    y[i, k, t] * (t + 1) <= C_ik[i, k],
                    name=f"Completion_j{i}_m{k}_t{t}",
                )

    # 6.-7. Sequencing: for each machine and each job pair, z_ij fixes the
    # order (z_ij=1 => i before j), via Big-M on the completion times.
    for i in range(n):
        for j in range(n):
            if i < j:
                for k in range(m):
                    model.addConstr(
                        C_ik[i, k] >= C_ik[j, k] - M * z_ij[i, j],
                        name=f"Seq1_j{i}_j{j}_m{k}",
                    )
                    model.addConstr(
                        C_ik[j, k] >= C_ik[i, k] - M * (1 - z_ij[i, j]),
                        name=f"Seq2_j{i}_j{j}_m{k}",
                    )

    # 8. Makespan.
    for i in range(n):
        model.addConstr(
            C_max >= C_ik[i, m - 1],
            name=f"Makespan_j{i}",
        )

    return model


def print_schedule(model: gp.Model, n: int, m: int, T: int) -> None:
    """Pretty-print the solved schedule, completion times and job order."""
    C_max = model.getVarByName("C_max")
    y = model.getVars()  # position-based lookup done via name below
    try:
        print(f"\n Makespan: {C_max.X:.2f}")
    except AttributeError:
        print("\n  Makespan: not available")

    var_by_name = {v.VarName: v for v in y}

    print("\n Job Schedule (y[i,k,t]=1):")
    for i in range(n):
        for k in range(m):
            active = []
            for t in range(T):
                v = var_by_name.get(f"y[{i},{k},{t}]")
                if v is not None and v.X > 0.5:
                    active.append(t + 1)
            if active:
                blocks = []
                start = end = active[0]
                for t in active[1:]:
                    if t == end + 1:
                        end = t
                    else:
                        blocks.append(f"{start}-{end}" if start != end else f"{start}")
                        start = end = t
                blocks.append(f"{start}-{end}" if start != end else f"{start}")
                print(f"  Job {i}, Machine {k}: Periods {', '.join(blocks)}")
            else:
                print(f"  Job {i}, Machine {k}: None")


def solve_instance(filename: str, time_limit: float = 900.0, mip_gap: float = 0.20) -> None:
    """Load one instance, build the model, solve it and print the schedule."""
    n, m, p_ik = read_instance(filename)
    print(f"\n{os.path.basename(filename)}: {n} jobs, {m} machines")

    model = build_flowshop_model(n, m, p_ik, time_limit, mip_gap)
    model.Params.TimeLimit = time_limit
    model.setParam("MIPGap", mip_gap)
    model.optimize()

    if model.status in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL):
        T = int(n * np.max(p_ik))
        print_schedule(model, n, m, T)
    else:
        print(f"\n Optimization status: {model.status}")  # e.g. INFEASIBLE


if __name__ == "__main__":
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    instances = sorted(glob.glob(os.path.join(data_dir, "VFR*.txt")))

    # Run all 10 VFR instances with a 15-minute time limit and 20% MIP gap
    # (weak enough for the 20-job instances to return feasible solutions).
    for instance in instances:
        solve_instance(instance, time_limit=900, mip_gap=0.20)