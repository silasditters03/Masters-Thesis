#!/usr/bin/env bash
# Run 3 MC, DCAFitterN, no PV refit.
# Full output goes to run.log; terminal shows reader progress.
#
# Before running:
# alienv enter O2/latest-834e08f85f-o2,JAliEn-ROOT/latest-834e08f85f-o2

set -euo pipefail

# ------------------------------------------------------------
# O2Physics development build
# ------------------------------------------------------------

ALICE_DIR="$HOME/SilasThesis/Masters-Thesis/alice"
O2PHYSICS_BUILD="$ALICE_DIR/sw/BUILD/O2Physics-latest-834e08f85f/O2Physics"

if [[ ! -d "$O2PHYSICS_BUILD/stage/bin" ]]; then
  printf 'ERROR: Build directory not found: %s\n' \
    "$O2PHYSICS_BUILD/stage/bin" >&2
  exit 1
fi

export PATH="$O2PHYSICS_BUILD/stage/bin:$PATH"

# ------------------------------------------------------------
# Input and configuration
# ------------------------------------------------------------

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

INPUT_FILE="${INPUT_FILE:-@$SCRIPT_DIR/input_files.txt}"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config_test.json}"
WORK_DIR="${WORK_DIR:-$SCRIPT_DIR/work-mc}"

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

# Only enable converters required by your AO2D table versions.
# Example:
# AOD_CONVERTERS='o2-analysis-mc-converter o2-analysis-mccollision-converter o2-analysis-tracks-extra-v002-converter'

converter_list="${AOD_CONVERTERS:-}"

if [[ -n "${converter_list//[[:space:]]/}" ]]; then
  read -r -a converters <<< "$converter_list"
  stages=("${converters[@]}" "${stages[@]}")
fi

# ------------------------------------------------------------
# Verify executables and inputs
# ------------------------------------------------------------

missing=0

for stage in "${stages[@]}"; do
  if ! command -v "$stage" >/dev/null 2>&1; then
    printf 'Missing executable: %s\n' "$stage" >&2
    missing=1
  fi
done

if (( missing )); then
  printf '%s\n' \
    'Load the compatible O2 runtime environment before running.' \
    'alienv enter O2/latest-834e08f85f-o2,JAliEn-ROOT/latest-834e08f85f-o2' >&2
  exit 127
fi

[[ -r "$INPUT_PATH" ]] || {
  printf 'Input is not readable: %s\n' "$INPUT_PATH" >&2
  exit 2
}

[[ -r "$CONFIG_FILE" ]] || {
  printf 'Configuration is not readable: %s\n' "$CONFIG_FILE" >&2
  exit 2
}

# ------------------------------------------------------------
# Working directory, full log and terminal progress
# ------------------------------------------------------------

mkdir -p -- "$WORK_DIR"
cd -- "$WORK_DIR"

exec > >(
  tee run.log | awk \
    -v input="$INPUT_PATH" \
    -v is_list="$input_is_list" \
    -v log_path="$PWD/run.log" '
    BEGIN {
      if (is_list) {
        while ((getline path < input) > 0) {
          sub(/\r$/, "", path)
          if (path !~ /^[[:space:]]*$/ && path !~ /^[[:space:]]*#/)
            wanted[path] = 1
        }
        close(input)
      } else {
        wanted[input] = 1
      }

      for (path in wanted) total++
      printf "AO2D reader: 0/%d files", total
      fflush()
    }

    /Read info:.*lfn=/ {
      path = $0
      sub(/^.*Read info:.*lfn=/, "", path)
      sub(/,.*/, "", path)

      if ((path in wanted) && !(path in seen)) {
        seen[path] = 1
        count++
        printf "\rAO2D reader: %d/%d files", count, total
        fflush()
      }
    }

    /^TCC_RUN_FINISHED=/ {
      split($0, status, "=")

      if (status[2] == 0)
        printf "\nWorkflow completed successfully.\n"
      else
        printf "\nWorkflow failed (exit code %s).\n", status[2]

      printf "Full log: %s\n", log_path
      fflush()
    }
  '
) 2>&1

args=(--configuration "json://$CONFIG_FILE" -b)

# These details are saved in run.log.
printf 'Input: %s\n' "$INPUT_FILE"
printf 'Configuration: %s\n' "$CONFIG_FILE"
printf 'Working directory: %s\n' "$PWD"
printf 'O2Physics build: %s\n' "$O2PHYSICS_BUILD"
printf 'Workflow: %s\n' "${stages[*]}"

# ------------------------------------------------------------
# Execute DPL workflow
# ------------------------------------------------------------

# The stdout pipes compose DPL workflow specifications.
# Only the final executable receives the AO2D reader options.

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

if run_stage 0 "$@"; then
  printf 'TCC_RUN_FINISHED=0\n'
else
  run_status=$?
  printf 'TCC_RUN_FINISHED=%s\n' "$run_status"
  exit "$run_status"
fi