#!/usr/bin/env bash
# Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
#
#
# Run the full PDK pipeline: build -> generate events -> fill TTree -> plot.
#
# Usage: ./run.sh [n_events] [seed] [model] [binding] [channel]
#   n_events : number of events to generate (default 10000)
#   seed     : RNG seed for reproducibility (default 1)
#   model    : nucleon momentum model (default polynomial)
#                polynomial|gfg|lfg|src|sf|benhar|ankowski
#   binding  : mean-field removal-energy model (default potential)
#                potential|constant|shell
#   channel  : decay channel (default pToKnu; LunarPDKGenerator --help lists all)
#
set -euo pipefail
cd "$(dirname "$0")"

N_EVENTS="${1:-10000}"
SEED="${2:-1}"
MODEL="${3:-polynomial}"
BINDING="${4:-potential}"
CHANNEL="${5:-pToKnu}"

echo ">> Building..."
make

echo ">> Generating $N_EVENTS '$CHANNEL' events ($MODEL, binding $BINDING, seed $SEED)..."
./build/LunarPDKGenerator --events "$N_EVENTS" --channel "$CHANNEL" \
   --model "$MODEL" --binding "$BINDING" --seed "$SEED" --fsi off \
   --config config/params.dat > data/output.txt

echo ">> Filling TTree (data/output.root)..."
root -l -b -q macros/text2tree.C

echo ">> Plotting distributions..."
root -l -b -q macros/plot_proton.C
root -l -b -q macros/plot_kaon.C
root -l -b -q macros/plot_removal.C

echo ">> Done. See plots/proton_p.png, plots/kaon_p.png, plots/removal_E.png"
