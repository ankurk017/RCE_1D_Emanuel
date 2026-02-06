# Emanuel RCE single-column model

1-D radiative–convective equilibrium model. Run with:

```bash
./run_RCE.sh --output_folder my_run --table params_ver2_table.in
```

## Using the shared environment

Activate the pre-built conda environment (no need to install conda or create an env yourself):

```bash
source /rtmp/akumar/conda/bin/activate rce
```

Then run the model and plotting as in the examples below. To recreate the env elsewhere, use `conda env create -f rce_environment.yml`.

## Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--output_folder` | Output directory (must not exist). | required |
| `--table` | Parameter table file (same format as `params_ver2_table.txt`). | `params_ver2_table.in` |
| `--plot yes` | After the run, generate plots into `output/figs/`. | no |

Positional form: you can pass the output folder as the first argument instead of `--output_folder NAME`.

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

Edit **`params_ver2_table.in`**  and **`sounding.in`** for SST. See **`params_ver2.in_README`** for parameter meanings; **`Users_guide.pdf`** for the full model description.

## Common issues

### "libgfortran.so.3: cannot open shared object file"

This can happen on a different machine (e.g. stratus) than the one the Linux executable was built on. Fix it by **rebuilding** the executable on this machine (requires `gfortran`):

1. **Get gfortran** (if needed). Try:
   - `module load gfortran`   (if available)
   - or `module spider gcc` / `module spider gfortran` and then e.g. `module load gcc/11.2.0`
   - or `module avail 2>&1 | grep -i gcc`

2. **Rebuild** the executable. Either use the Makefile:
   ```bash
   make rebuild
   ```
   Or build directly:
   ```bash
   gfortran -o rc_ver2unix rc_ver2.f
   ```


The new binary will run on this machine. If the build throws errors, **Ryan Driver or David Corredor** are the best people to help with library/compiler issues. If the gfortran error persists, contact **it@nsstc.uah.edu**.
