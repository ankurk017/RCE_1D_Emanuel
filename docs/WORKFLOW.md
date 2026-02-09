# Workflow: how to run the codebase

This diagram shows the **workflow the code follows** when you run it. As a user, you do **not** have to run all these steps yourself — you only need to activate the conda env and run `./run_RCE.sh` (and optionally edit input files). The flowchart is for **illustration and debugging** so you can see what happens under the hood.

```mermaid
flowchart TD
    A[Clone repo] --> B[Activate conda env: source .../activate rce]
    B --> C{Edit parameters?}
    C -->|Yes| D[Edit params_ver2_table.in and/or sounding.in]
    C -->|No| E[Run: ./run_RCE.sh --output_folder RUN_DIR --table FILE]
    D --> E
    E --> F[table_to_params.py: table → params_ver2.in]
    F --> G["execute_rce.sh: create run dir, symlink output"]
    G --> H["Run Fortran: rc_ver2unix"]
    H --> I{--plot / -p?}
    I -->|Yes| J[plot_rce_outputs.py → output/figs/]
    I -->|No| K[Outputs in run directory]
    J --> K
    K --> L[Done]
```

The flow includes: **execute_rce.sh** (create run dir, symlink `output`) → **Run Fortran** (`rc_ver2unix`).
