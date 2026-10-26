#!/usr/bin/env bash
# Batch job for one case on Tuolumne (LLNL, MI300A, flux). Submit with e.g.
#   flux batch -q pbatch -N 1 -B ufl -t 3h -n 4 -g 1 \
#     --output=build/<name>.out scripts/tuolumne_batch.sh <case dir> <casename> [nranks]
# The ppiclF library must already be built on a login node (gslib download):
#   FC=ftn CC=cc scripts/build_ppiclf.sh <case dir>/user_routines build/<case dir basename>
set -euo pipefail
case_dir=$1; casename=$2; nranks=${3:-4}
root=$(cd "$(dirname "$0")/.." && pwd)

module reset >/dev/null 2>&1 || true
module load PrgEnv-gnu gcc-native craype-accel-amd-gfx942 rocm/6.3.1 cmake python/3.12.2
module unload cray-libsci
export NEKRS_HOME=${NEKRS_HOME:-/p/lustre5/kim149/tuolumne/src/UF/nekrs-install}
export FC=ftn CC=cc
export LAUNCHER="flux run -g 1 -n"
cd "$root"
echo "== $(date)  $case_dir $casename on $nranks GPUs, $(hostname)"
scripts/run_case.sh "$case_dir" "$casename" "$nranks" --backend HIP --device-id 0
echo "== $(date)  done"
