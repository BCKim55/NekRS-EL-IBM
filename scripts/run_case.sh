#!/usr/bin/env bash
# Build the case's ppiclF library and run the case with nekRS.
#
#   scripts/run_case.sh <case dir> <casename> [nranks] [extra nekrs args]
#
# NEKRS_HOME must point at a nekRS install built from external/nekRS.
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
case_dir=$(cd "$1" && pwd)
casename=$2
nranks=${3:-1}
shift $(( $# < 3 ? $# : 3 ))

: "${NEKRS_HOME:?set NEKRS_HOME to the nekRS install}"
build="$root/build/$(basename "$case_dir")"
mkdir -p "$root/build"

if [ -f "$case_dir/user_routines/PPICLF_USER.h" ]; then
  "$root/scripts/build_ppiclf.sh" "$case_dir/user_routines" "$build" > "$build.ppiclf.log" 2>&1 \
    || { tail -30 "$build.ppiclf.log"; exit 1; }
  source "$build/udf.env"
else
  # device-only case (src/elgpu.hpp): no ppiclF
  export NEKRS_UDF_INCLUDES="$root/src"
  export NEKRS_UDF_LIBS=""
fi

cd "$case_dir"
# nekRS rebuilds the UDF only when the .udf changes, not when libppiclF.a does
rm -rf .cache/udf
mpirun -np "$nranks" "$NEKRS_HOME/bin/nekrs" --setup "$casename.par" "$@"
