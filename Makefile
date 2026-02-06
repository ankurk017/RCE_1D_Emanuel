# Build rc_ver2unix (Linux) from rc_ver2.f
# Use if the pre-built binary fails with "libgfortran.so.3: cannot open shared object file"
# (e.g. on a different machine like stratus). Requires gfortran: module load gcc, or install gcc-gfortran.

FC = gfortran
FFLAGS = -O2
TARGET = rc_ver2unix
# Static link Fortran runtime so the executable works without matching libgfortran on the system
STATIC = -static-libgfortran

$(TARGET): rc_ver2.f
	$(FC) $(FFLAGS) $(STATIC) -o $(TARGET) rc_ver2.f

clean:
	rm -f $(TARGET) *.o *.mod

# Force a full rebuild (use this when the existing binary was built on another machine)
rebuild: clean $(TARGET)

.PHONY: clean rebuild
