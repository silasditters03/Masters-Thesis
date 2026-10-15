#!/bin/bash
set -e
set -o pipefail

LOG_FILE="run.log"
exec > >(tee "$LOG_FILE") 2>&1

INPUT_FILE="data/001/AO2D.root"
ARGS="--configuration json://config.json -b"

echo "Launching native D0 MC workflow..."
o2-analysis-event-selection-service $ARGS | \
o2-analysis-propagationservice $ARGS | \
o2-analysis-trackselection $ARGS | \
o2-analysis-track-to-collision-associator $ARGS | \
o2-analysis-multcenttable $ARGS | \
o2-analysis-pid-tof-base $ARGS | \
o2-analysis-pid-tof $ARGS | \
o2-analysis-pid-tof-full $ARGS | \
o2-analysis-pid-tpc-service $ARGS | \
o2-analysis-hf-pid-creator $ARGS | \
o2-analysis-hf-track-index-skim-creator $ARGS | \
o2-analysis-hf-candidate-creator-2prong $ARGS | \
o2-analysis-hf-candidate-selector-d0 $ARGS | \
o2-analysis-hf-task-d0 $ARGS --aod-file $INPUT_FILE