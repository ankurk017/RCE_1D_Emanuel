# Emanuel RCE single-column model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue?style=flat-square)](https://www.python.org/downloads/)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey?style=flat-square)]()
[![Releases](https://img.shields.io/github/v/release/ankurk017/RCE_1D_Emanuel?include_prereleases&style=flat-square)](https://github.com/ankurk017/RCE_1D_Emanuel/releases)

1-D radiative–convective equilibrium model (Fortran core + Python for parameters and plotting).

This repository is an adaptation of the original **RCE MIT Single-Column Model** by **K. Emanuel** (15 March 2010). The Fortran core (`rc_ver2.f`) and scientific behavior are unchanged; the modifications are tooling and workflow around it:

- **Table-format parameters**: A human-readable parameter table (`params_ver2_table.in`) with Parameter, Value, Units, and Description columns. A Python script (`scripts/table_to_params.py`) converts it to the value-only `params_ver2.in` expected by the Fortran program.
- **Run and execute scripts**: `run_RCE.sh` provides a single entry point with options for output directory, parameter table file, and optional plotting. `execute_rce.sh` creates a user-specified run directory and a symlink so the model still reads/writes through `output/`, and picks the right executable by platform (Linux/Mac/Windows).
- **Python plotting**: `plotting/plot_rce_outputs.py` generates figures from model output (e.g. profiles, time series). Use `--plot yes` in `run_RCE.sh` to run the model and then produce plots in `output/figs/`.
- **Conda environment**: `rce_environment.yml` and a shared conda env so students can run without installing conda or building the Fortran executable themselves.
- **Build support**: A `Makefile` (and `make rebuild`) to recompile `rc_ver2unix` from `rc_ver2.f` on different systems (e.g. when `libgfortran` is missing).
- **Run-directory bookkeeping**: The run scripts copy the parameter table and a readable parameter summary into each run directory for reproducibility.
- **MATLAB scripts**: Parameter menus and run wrappers live in the `matlab/` folder (e.g. `run_model.m`, `get_parameters_ver2.m`). Run from the repo root with `addpath('matlab')` or from inside `matlab/`.

### Quick start

```bash
./run_RCE.sh --output_folder my_run --table params_ver2_table.in
```

## Using the shared environment

> [!NOTE]
> The scripts expect Python 3 with numpy and matplotlib. Activate the shared env (or recreate it from `rce_environment.yml`) before running.

```bash
source /rtmp/akumar/conda/bin/activate rce
```

Then run the model and plotting as in the examples below. To recreate the env elsewhere, use `conda env create -f rce_environment.yml`.

## Working on Meteor

To run the model on **Meteor**, connect via SSH (either directly or via a login node), then set up and build in the repo:

**Connect:**  
`ssh meteor`  
or, if required: `ssh stratus` or `ssh rayleigh`, then `ssh meteor`.

**Setup and test:**

```bash
# From the repo root after cloning
module load hdf5/1.14.6-gcc-11.5.0-cmrr
make rebuild
./run_RCE.sh --output_folder temp
```

> [!NOTE]
> The output directory (e.g. `temp`) must not already exist; the script will create it.

## Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--output_folder` | Output directory (must not exist). | required |
| `--table` | Parameter table file (same format as `params_ver2_table.txt`). | `params_ver2_table.in` |
| `--plot yes` | After the run, generate plots into `output/figs/`. | no |

Positional form: you can pass the output folder as the first argument instead of `--output_folder NAME`.

> [!NOTE]
> The output folder must not already exist; the script will create it. Choose a new name for each run (e.g. `SST25_run1`).

## Examples

```bash
# Default table, output in my_run/
./run_RCE.sh --output_folder my_run --table params_ver2_table.in

# Custom table file
./run_RCE.sh --output_folder SST28_run1 --table my_params.in

# Run and then generate plots
./run_RCE.sh --output_folder SST25_wind5 --table params_ver2_table.in --plot yes

# Positional output folder
./run_RCE.sh my_run --table params_ver2_table.in
./run_RCE.sh my_run
``` 

> [!NOTE]
> **Input files:** `params_ver2_table.in` is the parameter input file — edit it to change model parameters (e.g. SST, time step, radiation options). Also edit **`sounding.in`** for SST. See **`params_ver2.in_README`** for parameter meanings; **`Users_guide.pdf`** for the full model description.

> [!TIP]
> To get plots automatically after a run, add `--plot yes`; figures are written to `output/figs/`.

- **[Common issues](docs/TROUBLESHOOTING.md)** — libgfortran and other troubleshooting.
- **[Contributing](docs/CONTRIBUTING.md)** — How to contribute to this project.
- **[Workflow flowchart](docs/WORKFLOW.md)** — Diagram of how to run the codebase.

## Release

Stable versions are published as [Releases](https://github.com/ankurk017/RCE_1D_Emanuel/releases). Each release includes a tag (e.g. `v1.0.0`) and optional release notes. To use a specific version, clone the repo and check out the tag, or download the source tarball from the release page.

## Contact

For questions about this code or adaptation, contact **Ankur Kumar** at **ankur017@gmail.com** or **ankur.kumar@uah.edu**.

## License

This project is licensed under the [MIT License](./LICENSE).
