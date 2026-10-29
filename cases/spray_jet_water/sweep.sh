#!/usr/bin/env bash
# spray_jet_water sweep: coflow scaled with J^(1/3) so the entrainment demand of the
# droplet-driven jet is met without recirculation. Usage: sweep.sh p|i  (p = droplet-driven, i = gas-only inlet reference)
cd "$(dirname "$0")/../.."
export NEKRS_HOME=~/.local/nekrs OMPI_CXX=g++-16
kind=$1
run() { # idx phi dt steps t_avg u_co subcycles
  local name=$kind$1
  rm -rf build/spray/$name; rsync -aL --exclude .cache cases/spray_jet_water/ build/spray/$name/
  sed -i '' "s/^phi_inj = 5e-6$/phi_inj = $2/; s/^dt = 2.5e-3$/dt = $3/; s/^numSteps = 10000$/numSteps = $4/; s/^t_avg = 15.0$/t_avg = $5/; s/^u_co = 0.05$/u_co = $6/" build/spray/$name/jet.par
  [ "$7" -gt 1 ] && sed -i '' "s/^timeStepper = tombo2$/timeStepper = tombo2\nadvectionSubCyclingSteps = $7/" build/spray/$name/jet.par
  [ "$kind" = i ] && sed -i '' 's/^inlet_jet = 0$/inlet_jet = 1/' build/spray/$name/jet.par
  scripts/run_case.sh build/spray/$name jet 4 --backend CPU > build/spray_$name.log 2>&1
  echo "$name exit=$?" >> build/spray_sweep_$kind.out
}
run 1 5e-6 2.5e-3 10000 15.0 0.05 1
run 2 5e-5 1.25e-3 10000 7.5 0.11 1
run 3 5e-4 2e-4 25000 3.0 0.25 2
echo "sweep $kind done" >> build/spray_sweep_$kind.out
