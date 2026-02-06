#!/usr/bin/env bash
#
# Run the RCE model using a table-format parameter file.
# Converts params_ver2_table.in -> params_ver2.in, then runs execute_rce.sh.
#
# Usage:
#   ./run_RCE.sh --output_folder RUN_DIR [--table FILE] [--plot yes]
#   ./run_RCE.sh RUN_DIR [--table FILE] [--plot yes]   (positional still supported)
#
# Options:
#   --plot yes   After the run, generate plots into output/figs/ (default: no).
#
# Examples:
#   ./run_RCE.sh --output_folder test_output --table params_ver2_table.in
#   ./run_RCE.sh --output_folder SST25_fixedSST_wind5 --plot yes
#   ./run_RCE.sh SST28_run1 --table my_params_table.in --plot yes
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TABLE_FILE="params_ver2_table.in"
RUN_DIR=""
PLOT="no"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output_folder|--output-folder)
      [[ $# -ge 2 ]] || { echo "ERROR: --output_folder requires a path"; exit 2; }
      RUN_DIR="$2"
      shift 2
      ;;
    --table)
      [[ $# -ge 2 ]] || { echo "ERROR: --table requires a path"; exit 2; }
      TABLE_FILE="$2"
      shift 2
      ;;
    --plot)
      [[ $# -ge 2 ]] || { echo "ERROR: --plot requires a value (e.g. yes)"; exit 2; }
      if [[ "$2" == "yes" || "$2" == "y" ]]; then
        PLOT="yes"
      fi
      shift 2
      ;;
    -*)
      echo "ERROR: Unknown option: $1"
      echo "Usage: ./run_RCE.sh --output_folder RUN_DIR [--table FILE] [--plot yes]"
      echo "   or: ./run_RCE.sh RUN_DIR [--table FILE] [--plot yes]"
      exit 2
      ;;
    *)
      if [[ -z "$RUN_DIR" ]]; then
        RUN_DIR="$1"
        shift
      else
        echo "ERROR: Unexpected argument: $1"
        exit 2
      fi
      ;;
  esac
done

if [[ -z "$RUN_DIR" ]]; then
  echo "Usage: ./run_RCE.sh --output_folder RUN_DIR [--table FILE] [--plot yes]"
  echo "   or: ./run_RCE.sh RUN_DIR [--table FILE] [--plot yes]"
  exit 2
fi

if [[ ! -f "$SCRIPT_DIR/$TABLE_FILE" ]]; then
  echo "ERROR: Table file not found: $TABLE_FILE"
  echo "       Create it (same format as params_ver2_table.txt) or copy from a previous run."
  exit 1
fi

cd "$SCRIPT_DIR"
python3 "$SCRIPT_DIR/scripts/table_to_params.py" "$TABLE_FILE" params_ver2.in
"$SCRIPT_DIR/execute_rce.sh" "$RUN_DIR"

if [[ "$PLOT" == "yes" ]]; then
  if [[ -f "$SCRIPT_DIR/plotting/plot_rce_outputs.py" ]]; then
    echo "Generating plots..."
    python3 "$SCRIPT_DIR/plotting/plot_rce_outputs.py" --dir output --save output/figs
  else
    echo "WARNING: plotting/plot_rce_outputs.py not found; skipping plots."
  fi
fi
