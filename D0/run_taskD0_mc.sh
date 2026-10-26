
    
  
#!/usr/bin/env bash
# Run 3 MC, DCAFitterN, no PV refit, no reconstructed-centrality selection.
# Load a compatible O2Physics environment before running this script.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
INPUT_FILE="${INPUT_FILE:-@$SCRIPT_DIR/input_files.txt}"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.json}"

INPUT_PATH="$INPUT_FILE"
input_is_list=0

if [[ "$INPUT_PATH" == @* ]]; then
  input_is_list=1
  INPUT_PATH="${INPUT_PATH#@}"
fi

case "$INPUT_PATH" in /*) ;; *) INPUT_PATH="$PWD/$INPUT_PATH" ;; esac
case "$CONFIG_FILE" in /*) ;; *) CONFIG_FILE="$PWD/$CONFIG_FILE" ;; esac

if (( input_is_list )); then
  INPUT_FILE="@$INPUT_PATH"
else
  INPUT_FILE="$INPUT_PATH"
fi

WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/work-mc}"

# Only add converters required by the TABLE VERSIONS in this particular AO2D.
# Example (not a universal requirement):
# AOD_CONVERTERS='o2-analysis-mc-converter o2-analysis-mccollision-converter o2-analysis-tracks-extra-v002-converter'

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
  o2-analysis-hf-task-d0
)

converter_list="${AOD_CONVERTERS:-}"
if [[ -n "${converter_list//[[:space:]]/}" ]]; then
  read -r -a converters <<< "$converter_list"
  stages=("${converters[@]}" "${stages[@]}")
fi

missing=0
for stage in "${stages[@]}"; do
  if ! command -v "$stage" >/dev/null 2>&1; then
    printf 'Missing executable: %s\n' "$stage" >&2
    missing=1
  fi
done
if (( missing )); then
  echo 'Load an O2Physics installation compatible with the checked-in task sources.' >&2
  exit 127
fi
[[ -r "$INPUT_PATH" ]] || { printf 'Input is not readable: %s\n' "$INPUT_PATH" >&2; exit 2; }
[[ -r "$CONFIG_FILE" ]] || { printf 'Configuration is not readable: %s\n' "$CONFIG_FILE" >&2; exit 2; }
mkdir -p -- "$WORK_DIR"
cd -- "$WORK_DIR"
exec > >(tee run.log) 2>&1

args=(--configuration "json://$CONFIG_FILE" -b)
printf 'Input: %s\nConfiguration: %s\nWorking directory: %s\n' "$INPUT_FILE" "$CONFIG_FILE" "$PWD"
printf 'Workflow: %s\n' "${stages[*]}"

# The stdout pipes compose DPL workflow specifications. The final executable
# launches the combined graph and supplies the AO2D reader options once.
run_stage() {
  local index="$1"
  if (( index == ${#stages[@]} - 1 )); then
    "${stages[index]}" "${args[@]}" --run --aod-file "$INPUT_FILE" "${@:2}"
  else
    "${stages[index]}" "${args[@]}" | run_stage "$((index + 1))" "${@:2}"
  fi
}
run_stage 0 "$@"

