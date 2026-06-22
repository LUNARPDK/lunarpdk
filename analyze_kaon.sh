#!/usr/bin/env bash
# Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
#
#
# Generate events for every proton-momentum model and analyze the resulting
# lab-frame kaon momentum distribution. The K+ momentum is the experimentally
# observable p -> K+ nu signature: a free proton at rest gives a monochromatic
# kaon at p_K = (M_p^2 - m_K^2)/(2 M_p) ~ 0.339 GeV/c, which nuclear Fermi motion
# and binding (off-shell W) smear and shift downward.
#
# Usage: ./analyze_kaon.sh [n_events] [seed] [binding] [channel]
#   n_events : events per model (default 300000)
#   seed     : RNG seed (default 1)
#   binding  : mean-field removal-energy model (default potential)
#                potential|constant|shell
#   channel  : decay channel (default pToKnu; d2 = hadron daughter analysed)
#
set -euo pipefail
cd "$(dirname "$0")"

N_EVENTS="${1:-300000}"
SEED="${2:-1}"
BINDING="${3:-potential}"
CHANNEL="${4:-pToKnu}"
MODELS=(gfg lfg src sf hosm br gauss cfg benhar ankowski)

echo ">> Building..."
make

for m in "${MODELS[@]}"; do
   echo ">> Generating $N_EVENTS '$m' events ($CHANNEL, binding $BINDING)..."
   ./build/LunarPDKGenerator --events "$N_EVENTS" --channel "$CHANNEL" \
      --model "$m" --binding "$BINDING" --seed "$SEED" --fsi off > "data/kaon_${m}.txt"
done

echo ">> Analyzing kaon momentum distributions..."
root -l -b -q macros/plot_kaon_models.C

echo ">> Done. See plots/kaon_models.png and report/kaon_summary.txt"
