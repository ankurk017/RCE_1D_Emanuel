# Common issues

## "libgfortran.so.3: cannot open shared object file"

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
