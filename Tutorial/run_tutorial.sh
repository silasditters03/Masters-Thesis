#!/bin/bash
# Build and run one tutorial executable.
#
# Usage:   ./run_tutorial.sh <cmake-target> [extra task options]
# Example: ./run_tutorial.sh tutorial-associated-example
#
# Override the input file for one run with:
#   AOD_FILE=/path/to/other/AO2D.root ./run_tutorial.sh tutorial-task

set -euo pipefail

TARGET="${1:?Usage: $0 <cmake-target> [extra task options]}"
shift

# Project root = the folder this script lives in (it must contain CMakeLists.txt)
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$PROJECT_DIR/build"

AOD_FILE="${AOD_FILE:-$HOME/SilasThesis/Masters-Thesis/Tutorial/data/AO2D.root}"
SHM_SIZE=4000000000

# Configure once, then build only the requested target
if [ ! -f "$BUILD_DIR/CMakeCache.txt" ]; then
    cmake -S "$PROJECT_DIR" -B "$BUILD_DIR"
fi
cmake --build "$BUILD_DIR" --target "$TARGET" -j4

# Run from the build directory, where the executable lives
cd "$BUILD_DIR"

o2-analysis-tracks-extra-v002-converter -b --shm-segment-size "$SHM_SIZE" | \
o2-analysis-event-selection-service     -b --shm-segment-size "$SHM_SIZE" | \
o2-analysis-propagationservice          -b --shm-segment-size "$SHM_SIZE" | \
"./$TARGET"                             -b --shm-segment-size "$SHM_SIZE" \
    --aod-file "$AOD_FILE" "$@"