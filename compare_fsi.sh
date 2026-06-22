#!/usr/bin/env bash
# Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
#
#
# Compare the lab-frame meson momentum spectrum with and without final-state
# interactions (FSI), for the two benchmark channels:
#   p -> K+ nu    (kaon: nearly transparent in argon)
#   p -> e+ pi0   (pion: strongly absorbed / charge-exchanged near the Delta)
#
# Each channel is generated twice (--fsi off and --fsi on) with the same seed,
# then macros/compare_fsi.C overlays the two spectra, normalized per generated
# decay so the FSI-on curve shows both the depletion (smaller area = mesons lost
# to absorption / charge exchange) and the softening (elastic energy loss).
#
# Usage: ./compare_fsi.sh [n_events] [seed] [model]
#   n_events : events per channel per setting (default 200000)
#   seed     : RNG seed (default 1)
#   model    : nucleon-momentum model (default lfg)
#
set -euo pipefail
cd "$(dirname "$0")"

N_EVENTS="${1:-200000}"
SEED="${2:-1}"
MODEL="${3:-lfg}"

GEN=./build/LunarPDKGenerator

echo ">> Building..."
make >/dev/null
mkdir -p data plots

gen() {  # gen <channel> <fsi on|off> <outfile>
   "$GEN" --events "$N_EVENTS" --channel "$1" --model "$MODEL" --seed "$SEED" \
          --fsi "$2" > "$3" 2>/dev/null
}

echo ">> Generating p -> K+ nu  (FSI off / on, $N_EVENTS events, model $MODEL)..."
gen pToKnu  off data/fsi_pToKnu_off.txt
gen pToKnu  on  data/fsi_pToKnu_on.txt

echo ">> Generating p -> e+ pi0 (FSI off / on)..."
gen pToEPi0 off data/fsi_pToEPi0_off.txt
gen pToEPi0 on  data/fsi_pToEPi0_on.txt

echo ">> Plotting FSI on vs off overlay..."
root -l -b -q macros/compare_fsi.C

echo ">> Done: plots/fsi_compare.png"
