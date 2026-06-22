#!/usr/bin/env bash
# Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
#
#
# Compare the post-FSI meson outcome with the formation zone OFF vs ON, for the
# two benchmark channels:
#   p -> K+ nu    (kaon: slow, nearly transparent -> small formation length, mild shift)
#   p -> e+ pi0   (pion: strongly absorbed near the Delta -> large escape gain)
#
# Both runs keep the FSI cascade on (--fsi on) and use the same seed; the only
# difference is --formation-zone. A produced meson then free-streams a formation
# length L_f = (p/m)*c*tau_f before FSI is enabled, raising the escape ("none")
# fraction. The headline observable is the per-event primary-meson outcome, taken
# straight from the "# event N: ... outcome=<...>" headers.
#
# Usage: ./compare_fz.sh [n_events] [seed] [model] [c_tau_fm]
#   n_events : events per channel per setting (default 200000)
#   seed     : RNG seed (default 1)
#   model    : nucleon-momentum model (default lfg)
#   c_tau_fm : formation time c*tau_f in fm (default 0.342, SKAT)
#
set -euo pipefail
cd "$(dirname "$0")"

N_EVENTS="${1:-200000}"
SEED="${2:-1}"
MODEL="${3:-lfg}"
CTAU="${4:-0.342}"

GEN=./build/LunarPDKGenerator

echo ">> Building..."
make >/dev/null
mkdir -p data

gen() {  # gen <channel> <off|on> <outfile>  (FSI always on; off/on = formation zone)
    local fz_args=()
    [ "$2" = "on" ] && fz_args=(--formation-zone --formation-time "$CTAU")
    "$GEN" --events "$N_EVENTS" --channel "$1" --model "$MODEL" --seed "$SEED" \
           --fsi on "${fz_args[@]}" > "$3" 2>/dev/null
}

# Print the primary-meson outcome fractions from the per-event headers.
summary() {  # summary <file>
    awk '/^# event/ {
            for (i = 1; i <= NF; i++) if ($i ~ /^outcome=/) { split($i, a, "="); c[a[2]]++; n++ }
         }
         END {
            printf "    "
            split("none elastic cex produced absorbed", order, " ")
            for (k = 1; k <= 5; k++) { o = order[k]; if (o in c) printf "%s=%.1f%% ", o, 100.0*c[o]/n }
            printf "\n"
         }' "$1"
}

for ch in pToKnu pToEPi0; do
    echo ">> $ch  (FSI on; formation zone off vs on, c*tau=$CTAU fm, $N_EVENTS events, model $MODEL)..."
    gen "$ch" off "data/fz_${ch}_off.txt"
    gen "$ch" on  "data/fz_${ch}_on.txt"
    echo "  FZ off:"; summary "data/fz_${ch}_off.txt"
    echo "  FZ on :"; summary "data/fz_${ch}_on.txt"
done

echo ">> Done: data/fz_{pToKnu,pToEPi0}_{off,on}.txt"
echo "   (the post-FSI meson spectra can be overlaid with macros/compare_fsi.C,"
echo "    pointing it at the fz_*.txt pair instead of fsi_*.txt.)"
