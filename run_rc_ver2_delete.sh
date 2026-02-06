#!/usr/bin/env bash
#
# Wrapper that lets users choose an output directory.
# Picks executable by platform: rc_ver2mac (Mac), rc_ver2unix (Linux), rc_ver2dos.exe (Windows).
#
# Strategy:
# - Create a user-specified run directory (e.g. SST25_wind5)
# - Create/replace ./output as a symlink pointing to that directory
# - Ensure the target exists so `output/error.out` can be created
#
# Usage:
#   ./run_rc_ver2.sh RUN_DIR
#
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_rc_ver2.sh RUN_DIR

Notes:
  - RUN_DIR must not already exist (safer for students).ls
  - Executable is chosen by platform (rc_ver2mac / rc_ver2unix / rc_ver2dos.exe).
  - This script creates a symlink named `output` pointing to RUN_DIR.
  - The model reads inputs from:   sounding.in, params_ver2.in
  - The model writes outputs to:   output/*.out  (i.e., RUN_DIR/*.out)

Examples:
  ./run_rc_ver2.sh SST25_fixedSST_wind5
  ./run_rc_ver2.sh SST28_fixedSST_wind5
EOF
}

# Pick executable by platform
case "$(uname -s)" in
  Darwin)  exe="./rc_ver2mac" ;;
  Linux)   exe="./rc_ver2unix" ;;
  MINGW*|MSYS*|CYGWIN*) exe="./rc_ver2dos.exe" ;;
  *)       exe="./rc_ver2unix" ;;
esac
run_dir=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "ERROR: Unknown option: $1"
      usage
      exit 2
      ;;
    *)
      if [[ -z "$run_dir" ]]; then
        run_dir="$1"
        shift
      else
        echo "ERROR: Unexpected extra argument: $1"
        usage
        exit 2
      fi
      ;;
  esac
done

if [[ -z "$run_dir" ]]; then
  echo "ERROR: RUN_DIR is required."
  usage
  exit 2
fi

if [[ -e "$run_dir" ]]; then
  echo "ERROR: RUN_DIR already exists: $run_dir"
  echo "       Choose a new name (e.g., SST25_run1)."
  exit 2
fi

if [[ ! -f "$exe" ]]; then
  echo "ERROR: Executable not found: $exe"
  exit 2
fi

if [[ ! -x "$exe" ]]; then
  echo "ERROR: Executable is not marked executable: $exe"
  echo "       Try: chmod +x \"$exe\""
  exit 2
fi

echo "Run directory: $run_dir"
echo "Executable:    $exe"

cmds=()
cmds+=("mkdir -p \"$run_dir\"")

# If ./output exists:
# - if it's a symlink, remove it (safe)
# - if it's a real directory, move it into backup/ with a timestamped name
if [[ -L "output" ]]; then
  cmds+=("rm -f output")
elif [[ -e "output" ]]; then
  cmds+=("mkdir -p backup")
  backup_name="output_$(date +%Y%m%d_%H%M%S)"
  cmds+=("mv output \"backup/$backup_name\"")
  echo "Note: existing ./output will be moved to: backup/$backup_name"
fi

cmds+=("ln -s \"$run_dir\" output")

echo "Will run:"
printf '  %s\n' "${cmds[@]}"
echo "  \"$exe\""

eval "${cmds[0]}"
if [[ ${#cmds[@]} -ge 2 ]]; then eval "${cmds[1]}"; fi
if [[ ${#cmds[@]} -ge 3 ]]; then eval "${cmds[2]}"; fi
if [[ ${#cmds[@]} -ge 4 ]]; then eval "${cmds[3]}"; fi

# Copy params used for this run into output folder + pretty-table description
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_ABS="$(readlink -f output)"
if [[ -f "$SCRIPT_DIR/scripts/write_params_described.py" ]] && [[ -f "$SCRIPT_DIR/params_ver2.in" ]]; then
  python3 "$SCRIPT_DIR/scripts/write_params_described.py" "$SCRIPT_DIR/params_ver2.in" "$OUT_ABS"
fi

# Clean *.out files if any exist (fresh run semantics)
rm -f "output/"*.out 2>/dev/null || true

"$exe"

echo "Done."
echo "Outputs are in: $(readlink -f output)"

