# Flow-shop makespan MILP (Bowman) with Gurobi

**Short description:** Solving the permutation flow-shop scheduling problem to
minimise makespan using a time-indexed MILP (Bowman-style model) in Gurobi,
benchmarked on the VFR instances.

## What this is

Coursework for *Optimization Implementation in Production and Logistics* (OVGU
Magdeburg), Assignment 2. A flow-shop has `n` jobs that all visit `m` machines in
the same order; in a *permutation* flow-shop the jobs keep that order on every
machine. The goal is a schedule that minimises the **makespan** (the completion
time of the last job on the last machine).

The model is a classic **time-indexed / DP-based (Bowman) formulation**:

- Binary variable `y[i, k, t]` — job `i` occupies machine `k` at time `t`
- Continuous `C_ik[i, k]` — completion time of job `i` on machine `k`
- Binary `z_ij[i, j]` — sequencing order between job pairs
- Objective: minimise `C_max`

Constraints enforce machine capacity (one job per machine per slot), processing
times, precedence between successive machines, contiguous processing blocks, and
completion-time definition. Sequencing between jobs is Big-M enforced.

> An 8-constraint model with clear prose explanation is in
> `solve_flowshop.py:build_flowshop_model`; the notebook
> `Nasiya Pervez - Groubi_Assign.ipynb` is the original assignment submission.

## Repository layout

```
.
├── solve_flowshop.py            # Clean refactor: read instance -> build -> solve
├── Nasiya Pervez - Groubi_Assign.ipynb   # Original assignment notebook
├── data/
│   ├── VFR10_5_1_Gap.txt … VFR10_5_5_Gap.txt   # 10 jobs,  5 machines
│   └── VFR20_5_1_Gap.txt … VFR20_5_5_Gap.txt   # 20 jobs,  5 machines
└── LICENSE
```

## Data format

The VFR instances are from the classic flow-shop benchmark literature (VFR =
Vallada–Ruiz–Framinan), as provided by the lecturer. First line is
`[n jobs] [m machines]`; each following line is one job as
`[machine idx] [time] [machine idx] [time] …` (machines are 0-based).
See `data/ReadFirst.txt` equivalent description in the README of the source
session.

**Attribution:** the benchmark files are lecturer-supplied course material — they
are included here for reproducibility but are **not** covered by this repo's
license; original authorship/copyright remains with the dataset's creators.
Replace/add your own instances in `data/` if you use this elsewhere.

## Requirements & setup

- Python 3.9+
- `gurobipy` with a valid Gurobi license (academic licenses are free:
  https://www.gurobi.com/academia/academic-program-and-licenses/)
- `numpy`

```bash
pip install gurobipy numpy
python solve_flowshop.py
```

## How it works

`main` loops over every `data/VFR*.txt` instance and, for each one:

1. reads the processing-time matrix,
2. builds the model (horizon `T = n * max(p_ik)`, used as Big-M),
3. solves with a 15-minute time limit and a 20% MIP gap (the 20-job instances
   need that looser gap to return a feasible solution),
4. prints the makespan and a readable schedule.

Run a single file if you prefer:

```python
from solve_flowshop import solve_instance
solve_instance("data/VFR10_5_1_Gap.txt", time_limit=60)
```

## Expected results

On the small `VFR10_5_*` instances an optimal (or near-optimal) makespan is
found quickly; the `VFR20_5_*` instances are harder and hit the time limit,
returning a feasible solution within the 20% gap. Exact makespan values depend
on the Gurobi version and machine — the notebook header shows the bounds used.

## License

MIT — see [LICENSE](LICENSE). The VFR data files in `data/` are attributed,
not re-licensed.