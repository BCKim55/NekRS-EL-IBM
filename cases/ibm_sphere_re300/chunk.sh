#!/usr/bin/env bash
# One restart chunk of the Re 300 run, inside a flux allocation on Tuolumne
# (pdebug: 2 nodes = 8 GPUs, 60 min -> 2500 steps at 1.1 s/step; use 2000 for
# chunk 1, which also JIT-compiles). 20 000 steps = 8 chunks.
#
#   cases/ibm_sphere_re300/chunk.sh <k> [steps]     # k = 1, 2, ...
#
# Chunk k restarts from chunk_<k-1>.fld (the last checkpoint of chunk k-1,
# copied under a fixed name since nekRS restarts its file counter), and the
# force history goes to build/re300_chunk_<k>.out (time continues across chunks).
# Prebuild ppiclF on a login node first:
#   FC=ftn CC=cc scripts/build_ppiclf.sh cases/ibm_sphere_re300/user_routines build/ibm_sphere_re300
set -euo pipefail
cd "$(dirname "$0")/../.."
k=$1; steps=${2:-2500}; dir=cases/ibm_sphere_re300
if [ "$k" -gt 1 ]; then
  prev=chunk_$((k - 1)).fld
  [ -f "$dir/$prev" ] || { echo "missing $dir/$prev: run chunk $((k - 1)) first"; exit 1; }
  grep -q '^#\?startFrom' "$dir/sphere.par" || sed -i "s/^stopAt = numSteps$/startFrom = \"$prev\"\nstopAt = numSteps/" "$dir/sphere.par"
  sed -i "s|^#\?startFrom = .*|startFrom = \"$prev\"|" "$dir/sphere.par"
else
  sed -i 's|^startFrom = |#startFrom = |' "$dir/sphere.par"
fi
sed -i "s/^numSteps = .*/numSteps = $steps/; s/^checkpointInterval = .*/checkpointInterval = $steps/" "$dir/sphere.par"
grep -n "startFrom\|numSteps\|checkpointInterval" "$dir/sphere.par"
rm -f "$dir"/sphere0.f0*
REPO_ROOT=$PWD scripts/tuolumne_batch.sh "$dir" sphere 8 2>&1 | tee "build/re300_chunk_$k.out"
newest=$(ls -t "$dir"/sphere0.f0* | head -1)
cp -p "$newest" "$dir/chunk_$k.fld"
echo "chunk $k done: $newest -> $dir/chunk_$k.fld; last line: $(grep 'sphere: step' "build/re300_chunk_$k.out" | tail -1)"
