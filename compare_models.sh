#!/usr/bin/env bash
# Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
#
#
# Generate events for every nuclear momentum model and overlay their nucleon
# momentum distributions.
#
# Usage: ./compare_models.sh [n_events] [seed] [binding] [channel]
#   n_events : events per model (default 200000)
#   seed     : RNG seed (default 1)
#   binding  : mean-field removal-energy model (default potential)
#                potential|constant|shell
#   channel  : decay channel (default pToKnu)
#
set -euo pipefail
cd "$(dirname "$0")"

N_EVENTS="${1:-200000}"
SEED="${2:-1}"
BINDING="${3:-potential}"
CHANNEL="${4:-pToKnu}"
MODELS=(gfg lfg src sf hosm br gauss cfg benhar ankowski)

echo ">> Building..."
make

for m in "${MODELS[@]}"; do
   echo ">> Generating $N_EVENTS '$m' events ($CHANNEL, binding $BINDING)..."
   ./build/LunarPDKGenerator --events "$N_EVENTS" --channel "$CHANNEL" \
      --model "$m" --binding "$BINDING" --seed "$SEED" --fsi off > "data/proton_${m}.txt"
done

echo ">> Overlaying proton momentum distributions..."
root -l -b -q macros/plot_models.C

echo ">> Done. See plots/proton_models.png"
