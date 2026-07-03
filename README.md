# Dengue Heterogeneous Agent-Based Model (HABM) — Santa Marta, Colombia

This repository contains the source code, input data, calibration/analysis
scripts, and manuscript for a spatially explicit **agent-based model of dengue
transmission** in Santa Marta, Colombia. The model couples a distributed C++
simulation (built on [Repast HPC](https://repast.github.io/)) with Python
scripts for calibration and sensitivity analysis.

The simulation represents **human movement** across an urban spatial grid and
**mosquito (*Aedes aegypti*) population dynamics** using a per-patch
susceptible–exposed–infectious (SEI) compartmental model, driven by
neighbourhood-level entomological (House Index) data.

> **Note for reviewers/editors.** This is a research *simulator*: it is designed
> to reproduce observed epidemiological patterns for Santa Marta, not to forecast
> future incidence. The compiled manuscript is in [`Draft_DENGUE/`](Draft_DENGUE/).

---

## Repository structure

```
.
├── README.md                     ← you are here
├── QUICKSTART.md                 ← condensed build-and-run instructions
├── SENSITIVITY_ANALYSIS_README.md← sensitivity-analysis methodology & usage
│
├── C++ simulation (Repast HPC)
│   ├── Main.cpp                  ← entry point; reads props, builds Params
│   ├── MyModel.{h,cpp}           ← top-level model: grid, layers, scheduling
│   ├── SEIModel.{h,cpp}          ← per-patch mosquito SEI dynamics
│   ├── MyHuman.{h,cpp}           ← human agent movement & infection
│   ├── MyReadData.{h,cpp}        ← loads location CSVs & temperature series
│   ├── CSVWriter.{h,cpp}         ← writes simulation output CSVs
│   ├── Makefile                  ← build rules
│   └── env                       ← compiler/library paths (edit for your system)
│
├── props/                        ← simulation configuration
│   ├── model.props               ← model parameters (population, grid, rates, HI)
│   ├── config.props              ← Repast HPC logging config
│   └── model_baseline.props      ← baseline parameter set
│
├── Python workflows
│   ├── calibration_I0.py         ← calibrate initial infected count (I0)
│   ├── calibration_bo.py         ← Bayesian-optimisation calibration
│   ├── run_sensitivity.py        ← one-at-a-time sensitivity analysis
│   ├── run_sa_corrected.py       ← corrected-model sensitivity runs
│   ├── test_params.py            ← parameter-parsing sanity check
│   ├── pyproject.toml / uv.lock  ← Python environment (managed with uv)
│   └── *.ipynb                   ← calibration & sensitivity notebooks
│
├── Input data
│   ├── new_grid.csv, santa_marta_grid_relative.csv ← 135×200 spatial grid
│   ├── locations.csv             ← residential/work/study/other locations
│   ├── hi_by_neighborhood_2023.csv, hi_activation_map.csv ← House Index data
│   ├── temperatura-precipitacion-bello.csv         ← temperature time series
│   ├── *.geojson                 ← neighbourhood / grid geometries
│   └── santa_marta_data/         ← raw epidemiological Excel inputs (2023–2024)
│
├── view/                         ← analysis notebooks, processed data & figures
│
└── Draft_DENGUE/                 ← the manuscript (see below)
    ├── template_new.tex          ← LaTeX source (MDPI / Applied Sciences)
    ├── template_new.pdf          ← compiled manuscript
    ├── references.bib            ← bibliography
    ├── Images/                   ← manuscript figures
    └── Definitions/              ← MDPI class & bibliography-style files
```

## The manuscript

The paper lives in [`Draft_DENGUE/`](Draft_DENGUE/) and is self-contained: the
MDPI document class (`Definitions/mdpi.cls`), bibliography styles, figures, and
`references.bib` are all included, so `template_new.tex` compiles as-is with a
standard LaTeX + BibTeX toolchain. `template_new.pdf` is the compiled version.

## Requirements

**C++ simulation**
- [Repast HPC 2.3.1](https://repast.github.io/)
- OpenMPI, Boost
- A C++17 compiler

Compiler and library paths are set in the [`env`](env) file — edit these to
point at your local Repast HPC / Boost / MPI installation before building.

**Python workflows** — Python ≥ 3.10. Dependencies are pinned in
[`pyproject.toml`](pyproject.toml) / [`uv.lock`](uv.lock) and are most easily
installed with [uv](https://github.com/astral-sh/uv):

```bash
uv sync
```

## Building & running the simulation

```bash
# Build
make clean && make

# Run a single-process simulation
mpirun -n 1 ./Main.exe ./props/config.props ./props/model.props
```

Key parameters (population size, grid, disease rates, House-Index activation)
are set in [`props/model.props`](props/model.props). See
[`QUICKSTART.md`](QUICKSTART.md) for full details on the parameter flow and
model architecture.

## Calibration

Find the best initial infected count `I0`:

```bash
python3 calibration_I0.py --mode run --r-start 0 --r-end 10 --replicates 10
python3 calibration_I0.py --mode aggregate
python3 calibration_I0.py --mode final --replicates 10
```

A Bayesian-optimisation variant (over `I0`, bite-aggregation `z`, and
transmission `beta_mh`) is provided in `calibration_bo.py`.

## Sensitivity analysis

One-at-a-time sensitivity analysis with Common Random Numbers:

```bash
python3 run_sensitivity.py --mode run --replicates 3 --outdir out_sa
python3 run_sensitivity.py --mode aggregate --outdir out_sa
```

See [`SENSITIVITY_ANALYSIS_README.md`](SENSITIVITY_ANALYSIS_README.md) for the
methodology (normalised sensitivity coefficients, perturbed parameters, etc.).

## Reproducibility notes

- Simulation run outputs (per-seed CSV/NumPy ensembles under `*_runs/`,
  `test_run/`, etc.) are **regenerable** from the scripts above and are therefore
  not committed; only scripts and key summary results are tracked.
- Random seeds are derived deterministically as `base_seed + replicate_id`, so
  using the same replicate ids across configurations synchronises the random
  streams for fair comparison.

## License / contact

Corresponding author: Sara Cano — *smcanom* on GitHub.
Please open an issue for questions about building or running the model.
