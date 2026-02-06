#!/usr/bin/env python3
"""
Write a copy of params_ver2.in into the output folder and a pretty-table
description file (params_ver2_table.txt) for that run.

Usage:
  python write_params_described.py PARAMS_FILE OUTPUT_DIR
  python write_params_described.py params_ver2.in output
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

# Order of value lines in params_ver2.in (1-based line numbers in file).
# Blank/header lines are skipped by rc_ver2.f; only these lines carry a value.
VALUE_LINE_NUMBERS = [
    4, 5, 6, 7, 8,   # integration
    11, 12, 13, 14, 15, 16, 17, 18, 19,  # radiation/calendar
    20, 21, 22, 23, 24, 25,
    28, 29, 30, 31, 32, 33, 34, 35,  # greenhouse gases + multipliers
    38, 39,   # convection
    42, 43, 44, 45,   # surface fluxes
    48, 49, 50, 51, 52, 53,   # omega
    56, 57,   # WTG
]

# (parameter name, short description, units) for each value line, in same order.
PARAM_META = [
    ("RESTART", "Restart from previous run? (y: use sounding.out; n: use sounding.in)", "y/n"),
    ("ENDTIME", "End time of integration", "days"),
    ("DT", "Time step", "minutes"),
    ("AVTIME", "Averaging time for final profiles", "days"),
    ("PRINFREQ", "Output sampling interval for time series", "hours"),
    ("RADINT", "Interactive radiation?", "y/n"),
    ("ICLDS", "Interactive clouds?", "y/n"),
    ("TSINT", "Interactive surface temperature (ocean)? n = fixed SST", "y/n"),
    ("RADFREQ", "Frequency of radiation calls", "hours"),
    ("SCON", "Solar constant", "W m^-2"),
    ("RLAT", "Latitude", "degrees"),
    ("MONTH", "Starting month", "1-12"),
    ("IDAY", "Starting day", "1-31"),
    ("HOUR", "Starting hour", "0-23"),
    ("TDEP", "Time-dependent radiation?", "y/n"),
    ("DDEP", "Date-dependent radiation?", "y/n"),
    ("DARAD", "Diurnal-average radiation?", "y/n"),
    ("ANRAD", "Annual-average radiation?", "y/n"),
    ("CALB", "Calculate ocean albedo?", "y/n"),
    ("ALB", "Fixed surface albedo (if CALB=n)", "0-1"),
    ("CO2", "CO2 concentration", "ppm"),
    ("CH4", "CH4 concentration", "ppm"),
    ("N2O", "N2O concentration", "ppb"),
    ("CFC11", "CFC-11 concentration", "ppt"),
    ("CFC12", "CFC-12 concentration", "ppt"),
    ("H2OI", "Interactive water vapor?", "y/n"),
    ("H2OM", "Radiation H2O multiplier", "-"),
    ("O3M", "Radiation O3 multiplier", "-"),
    ("DCONV", "Dry adiabatic adjustment?", "y/n"),
    ("MCONV", "Moist convection scheme?", "y/n"),
    ("FLUXSWITCH", "Turbulent surface fluxes on/off", "y/n"),
    ("BETA", "Fraction of surface covered by water", "0-1"),
    ("DEPTH", "Mixed-layer depth", "m"),
    ("VS0", "Surface wind speed", "m s^-1"),
    ("WCUBE", "Use cubic profile of omega?", "y/n"),
    ("WMAX", "Extreme value of omega (negative = upward)", "hPa hour^-1"),
    ("PERIOD", "Period of omega variation", "days"),
    ("PWZERO", "Pressure where omega = 0 (top)", "hPa"),
    ("PWBOTTOM", "Pressure where omega = 0 (bottom)", "hPa"),
    ("PWMAX", "Pressure where omega reaches extreme", "hPa"),
    ("WTG", "Apply weak-temperature-gradient approximation?", "y/n"),
    ("PFIX", "Pressure above which sounding fixed", "hPa"),
]


def main() -> None:
    if len(sys.argv) != 3:
        print("Usage: write_params_described.py PARAMS_FILE OUTPUT_DIR", file=sys.stderr)
        sys.exit(2)
    params_file = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])

    if not params_file.is_file():
        print(f"ERROR: Not a file: {params_file}", file=sys.stderr)
        sys.exit(1)
    out_dir.mkdir(parents=True, exist_ok=True)

    lines = params_file.read_text().splitlines()
    # 1-based line numbers -> 0-based index
    value_indices = [n - 1 for n in VALUE_LINE_NUMBERS]
    if len(PARAM_META) != len(value_indices):
        print("ERROR: PARAM_META length does not match VALUE_LINE_NUMBERS", file=sys.stderr)
        sys.exit(1)

    values = []
    for i in value_indices:
        if i < len(lines):
            values.append(lines[i].strip())
        else:
            values.append("")

    if len(values) != len(PARAM_META):
        print("ERROR: extracted value count mismatch", file=sys.stderr)
        sys.exit(1)

    # 1) Exact copy into output folder
    dest_params = out_dir / "params_ver2.in"
    shutil.copy2(params_file, dest_params)

    # 2) Pretty table
    dest_table = out_dir / "params_ver2_table.txt"
    col_name = "Parameter"
    col_value = "Value"
    col_units = "Units"
    col_desc = "Description"
    w_name = max(len(col_name), max(len(m[0]) for m in PARAM_META))
    w_value = max(len(col_value), max(len(v) for v in values))
    w_units = max(len(col_units), max(len(m[2]) for m in PARAM_META))
    w_desc = max(len(col_desc), max(len(m[1]) for m in PARAM_META))

    sep = "  "
    header = f"{col_name:<{w_name}}{sep}{col_value:<{w_value}}{sep}{col_units:<{w_units}}{sep}{col_desc}"
    ruler = "-" * len(header)

    rows = [header, ruler]
    for (name, desc, units), val in zip(PARAM_META, values):
        rows.append(f"{name:<{w_name}}{sep}{val:<{w_value}}{sep}{units:<{w_units}}{sep}{desc}")

    dest_table.write_text("\n".join(rows) + "\n")
    print(f"Wrote {dest_params} and {dest_table}")


if __name__ == "__main__":
    main()
