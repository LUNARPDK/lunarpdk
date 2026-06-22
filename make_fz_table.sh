#!/usr/bin/env bash
# Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
#
#
# Tabulate the formation-zone impact on the primary-meson survival yield, for the
# six benchmark FSI channels, as a sensitivity scan over the formation time
# c*tau_f. The baseline (no formation zone) reuses the FSI-on samples already
# produced by make_plots.sh (data/fsi_<ch>_on.txt); the +FZ columns are generated
# here with --formation-zone --formation-time <c*tau_f>. Everything uses the same
# statistics as Table tab:fsi in the report (lfg, seed 1, 2e5 events) so the
# baseline column reproduces the tabulated FSI-on yields.
#
# Yield = survival probability = fraction of decays whose primary meson ends as
# none + elastic + produced (the epsilon_FSI used in the event-rate prediction);
# "none" is the clean-escape (untouched) fraction.
#
# Usage: ./make_fz_table.sh [n_events] [seed] [model]
#   n_events : events per channel per setting (default 200000)
#   seed     : RNG seed (default 1)
#   model    : nucleon-momentum model (default lfg)
#
# Writes report/fz_summary.txt (sourced by tab:fz in report.tex / lunar_paper.tex).
set -euo pipefail
cd "$(dirname "$0")"

N="${1:-200000}"
SEED="${2:-1}"
MODEL="${3:-lfg}"
CTAUS=(0.342 0.5 1.0)
GEN=./build/LunarPDKGenerator

# Channels in tab:fsi order, with their meson label.
CHANS=(pToKnu pToMuK0 pToEEta pToEPi0 pToNuPip nToEKm)
MESONS=("K+" "K0" "eta" "pi0" "pi+" "K-")

echo ">> Building..."
make >/dev/null
mkdir -p data report

# yield_none <file> -> "<yield> <none>" as fractions of all classified events.
yield_none() {
    awk '/^# event/ {
            p = index($0, "outcome="); if (!p) next
            oc = substr($0, p + 8); sub(/[ \t\r\n]+$/, "", oc)
            n++
            if (oc == "none" || oc == "elastic" || oc == "produced") y++
            if (oc == "none") none++
         }
         END { if (n) printf "%.4f %.4f", y/n, none/n; else printf "NA NA" }' "$1"
}

OUT=report/fz_summary.txt
{
    echo "# Formation-zone impact on primary-meson survival yield (and clean-escape 'none' fraction)."
    echo "# Same sample as tab:fsi: model=$MODEL seed=$SEED events=$N. Baseline = FSI on, no formation zone."
    echo "# c*tau_f columns in fm: ${CTAUS[*]}"
    printf "# %-9s %6s | %-30s | %-30s\n" "channel" "meson" "yield: base ${CTAUS[*]}" "none:  base ${CTAUS[*]}"
} > "$OUT"

for i in "${!CHANS[@]}"; do
    ch="${CHANS[$i]}"; meson="${MESONS[$i]}"
    base="data/fsi_${ch}_on.txt"
    if [ ! -f "$base" ]; then
        echo ">> [$ch] baseline missing; generating..."
        "$GEN" --events "$N" --channel "$ch" --model "$MODEL" --seed "$SEED" --fsi on \
               > "$base" 2>/dev/null
    fi

    echo ">> [$ch] ($meson) tallying baseline + formation zone scan..."
    read -r y_base n_base <<<"$(yield_none "$base")"

    ys=("$y_base"); ns=("$n_base")
    for ct in "${CTAUS[@]}"; do
        f="data/fz_${ch}_ct${ct}.txt"
        "$GEN" --events "$N" --channel "$ch" --model "$MODEL" --seed "$SEED" \
               --fsi on --formation-zone --formation-time "$ct" > "$f" 2>/dev/null
        read -r y n <<<"$(yield_none "$f")"
        ys+=("$y"); ns+=("$n")
    done

    printf "  %-9s %6s | %s | %s\n" "$ch" "$meson" "${ys[*]}" "${ns[*]}" >> "$OUT"
done

echo ">> Done. Summary:"
cat "$OUT"
