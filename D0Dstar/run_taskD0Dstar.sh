
#!/usr/bin/env bash
# Run 3 MC, DCAFitterN, no PV refit, no reconstructed-centrality selection.
# Uses the existing O2Physics 834e08f85f development build.
#
# Before running, enter:
#   alienv enter ninja/latest-834e08f85f-o2

set -euo pipefail

# ------------------------------------------------------------
# O2Physics development build setup
# ------------------------------------------------------------

ALICE_DIR="$HOME/SilasThesis/Masters-Thesis/alice"

O2PHYSICS_BUILD="$ALICE_DIR/sw/BUILD/O2Physics-latest-834e08f85f/O2Physics"

# Check that the compiled O2Physics executables exist.
if [[ ! -d "$O2PHYSICS_BUILD/stage/bin" ]]; then
  echo "ERROR: O2Physics build directory not found:"
  echo "$O2PHYSICS_BUILD/stage/bin"
  exit 1
fi

# Make the executables available to this script.
export PATH="$O2PHYSICS_BUILD/stage/bin:$PATH"

echo "Using O2Physics build:"
echo "$O2PHYSICS_BUILD"

# ------------------------------------------------------------
# Input and configuration
# ------------------------------------------------------------

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
INPUT_FILE="${INPUT_FILE:-@$SCRIPT_DIR/input_files.txt}"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config_test.json}"

INPUT_PATH="$INPUT_FILE"
input_is_list=0

if [[ "$INPUT_PATH" == @* ]]; then
  input_is_list=1
  INPUT_PATH="${INPUT_PATH#@}"
fi

case "$INPUT_PATH" in
  /*) ;;
  *) INPUT_PATH="$PWD/$INPUT_PATH" ;;
esac

case "$CONFIG_FILE" in
  /*) ;;
  *) CONFIG_FILE="$PWD/$CONFIG_FILE" ;;
esac

if (( input_is_list )); then
  INPUT_FILE="@$INPUT_PATH"
else
  INPUT_FILE="$INPUT_PATH"
fi

WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/work-mc}"

# ------------------------------------------------------------
# Optional AOD converters
# ------------------------------------------------------------

# Only add converters required by the TABLE VERSIONS
# in this particular AO2D.
#
# Example:
# AOD_CONVERTERS='o2-analysis-mc-converter o2-analysis-mccollision-converter o2-analysis-tracks-extra-v002-converter'

# ------------------------------------------------------------
# Analysis pipeline
# ------------------------------------------------------------

stages=(
  o2-analysis-event-selection-service
  o2-analysis-propagationservice
  o2-analysis-trackselection
  o2-analysis-track-to-collision-associator
  o2-analysis-multcenttable
  o2-analysis-lf-mc-centrality
  o2-analysis-ft0-corrected-table
  o2-analysis-pid-tof-merge
  o2-analysis-pid-tpc-service
  o2-analysis-hf-pid-creator
  o2-analysis-hf-track-index-skim-creator
  o2-analysis-hf-candidate-creator-2prong
  o2-analysis-hf-candidate-selector-d0
  o2-analysis-hf-candidate-creator-dstar
  o2-analysis-hf-candidate-selector-dstar-to-d0-pi
  o2-analysis-hf-task-d0-dstar
)

converter_list="${AOD_CONVERTERS:-}"

if [[ -n "${converter_list//[[:space:]]/}" ]]; then
  read -r -a converters <<< "$converter_list"
  stages=("${converters[@]}" "${stages[@]}")
fi

# ------------------------------------------------------------
# Verify executables
# ------------------------------------------------------------

missing=0

for stage in "${stages[@]}"; do
  if ! command -v "$stage" >/dev/null 2>&1; then
    printf 'Missing executable: %s\n' "$stage" >&2
    missing=1
  fi
done

if (( missing )); then
  echo "ERROR: Some O2Physics executables are missing." >&2
  echo "Enter the matching Ninja environment:" >&2
  echo "alienv enter ninja/latest-834e08f85f-o2" >&2
  exit 127
fi

# ------------------------------------------------------------
# Verify input files
# ------------------------------------------------------------

[[ -r "$INPUT_PATH" ]] || {
  printf 'Input is not readable: %s\n' "$INPUT_PATH" >&2
  exit 2
}

[[ -r "$CONFIG_FILE" ]] || {
  printf 'Configuration is not readable: %s\n' "$CONFIG_FILE" >&2
  exit 2
}

# ------------------------------------------------------------
# Working directory and logging
# ------------------------------------------------------------

mkdir -p -- "$WORK_DIR"
cd -- "$WORK_DIR"

exec > >(tee run.log) 2>&1

args=(--configuration "json://$CONFIG_FILE" -b)

printf 'Input: %s\n' "$INPUT_FILE"
printf 'Configuration: %s\n' "$CONFIG_FILE"
printf 'Working directory: %s\n' "$PWD"
printf 'O2Physics build: %s\n' "$O2PHYSICS_BUILD"
printf 'Workflow: %s\n' "${stages[*]}"

# ------------------------------------------------------------
# Execute DPL workflow
# ------------------------------------------------------------

# The stdout pipes compose DPL workflow specifications.
# The final executable launches the combined graph and
# supplies the AO2D reader options once.

run_stage() {
  local index="$1"

  if (( index == ${#stages[@]} - 1 )); then

    "${stages[index]}" "${args[@]}" --run --aod-file "$INPUT_FILE" \
      --shm-segment-size 4000000000 \
      --aod-memory-rate-limit 1000000000 \
      "${@:2}"

  else

    "${stages[index]}" "${args[@]}" |
      run_stage "$((index + 1))" "${@:2}"

  fi
}

run_stage 0 "$@"
