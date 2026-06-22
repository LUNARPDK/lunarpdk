#!/usr/bin/env bash
# Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
#
#
# Rebuild every figure in plots/ from scratch with the shared styling
# (macros/pdk_style.h). Regenerates all the data each plot needs, then runs each
# analysis macro. One command to refresh the whole figure set consistently.
#
# Usage: ./make_plots.sh [n_events] [seed] [channel]
#   n_events : events per sample (default 200000)
#   seed     : RNG seed (default 1)
#   channel  : decay channel (default pToKnu; LunarPDKGenerator --help lists all)
#
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p plots report data

N="${1:-200000}"
SEED="${2:-1}"
CHANNEL="${3:-pToKnu}"
GEN=./build/LunarPDKGenerator

# Wrapper: gen <model> <binding> -> stdout, for the active channel/seed/N.
gen() { "$GEN" --events "$N" --channel "$CHANNEL" --model "$1" --binding "$2" \
              --seed "$SEED" --fsi off --config config/params.dat; }

echo ">> Building..."
make

# 1) Single-distribution plots (nucleon_p, d2_p, removal_E) from a TTree.
echo ">> [1/6] single-model sample -> TTree"
gen lfg potential > data/output.txt
root -l -b -q macros/text2tree.C

# 2) Nucleon momentum overlay across all models.
echo ">> [2/6] nucleon spectra across models"
for m in gfg lfg src sf hosm br gauss cfg benhar ankowski; do
   gen "$m" potential > "data/proton_${m}.txt"
done

# 3) Hadron-daughter momentum overlay across all models.
echo ">> [3/6] daughter spectra across models"
for m in gfg lfg src sf hosm br gauss cfg benhar ankowski; do
   gen "$m" potential > "data/kaon_${m}.txt"
done

# 4) Tabulated-SF removal energy: argon proton grid vs neutron grid.
echo ">> [4/6] tabulated-SF removal energy (proton vs neutron grid)"
"$GEN" --events "$N" --channel pToKnu  --model benhar --seed "$SEED" --fsi off > data/erem_sf_proton.txt
"$GEN" --events "$N" --channel nToNuK0 --model benhar --seed "$SEED" --fsi off > data/erem_sf_neutron.txt

# 5) Daughter spectra across binding models (two analytic momentum models;
#    binding does not apply to the tabulated benhar/ankowski models).
echo ">> [5/6] daughter spectra across binding models"
for m in lfg sf; do for b in potential constant shell; do
   gen "$m" "$b" > "data/kbind_${m}_${b}.txt"
done; done

# All 14 decay channels in the registry (8 proton + 6 neutron modes).
ALL_CHANNELS="pToKnu pToEPi0 pToMuPi0 pToNuPip pToEEta pToMuEta pToEK0 pToMuK0 \
              nToEPim nToMuPim nToNuPi0 nToNuEta nToNuK0 nToEKm"
# The six benchmark channels that get a full-statistics FSI on-vs-off spectrum
# figure (K+ transparent; pi+/pi0 strongly interacting; eta; K0; K- antikaon).
FSI_SPECTRUM="pToKnu pToNuPip pToEPi0 pToEEta pToMuK0 nToEKm"
# Reduced count for the per-channel/per-model efficiency sweeps below: these feed
# fractions (window and FSI yields), so 50k events is plenty and keeps runtime sane.
NE=50000

# 6) FSI study: full-statistics spectrum overlays for the six benchmark channels
#    (off + on), plus an FSI-on sample for every channel so fsi_outcomes.C can
#    tabulate the primary-meson fate / survival yield for all 14 modes.
echo ">> [6/8] FSI spectra (6 benchmark channels, off+on @ $N)"
for ch in $FSI_SPECTRUM; do
   for f in off on; do
      "$GEN" --events "$N" --channel "$ch" --model lfg --seed "$SEED" --fsi "$f" \
             > "data/fsi_${ch}_${f}.txt"
   done
done
echo ">> [7/8] FSI-on yields for the remaining channels (@ $NE)"
for ch in $ALL_CHANNELS; do
   [ -f "data/fsi_${ch}_on.txt" ] && continue   # already produced above
   "$GEN" --events "$NE" --channel "$ch" --model lfg --seed "$SEED" --fsi on \
          > "data/fsi_${ch}_on.txt"
done

# 8) Window-efficiency sweep: hadron-daughter momentum for every channel and
#    every nuclear model (FSI off, legacy table), feeding window_eff.C. The
#    generator's structured "# forbidden ..." stderr summary is collected here
#    into report/forbidden_log.txt for offshell_W.C (the forbidden-fraction bar).
echo ">> [8/8] window-efficiency sweep (14 channels x 10 models @ $NE)"
: > report/forbidden_log.txt
for ch in $ALL_CHANNELS; do
   for m in gfg lfg src sf hosm br gauss cfg benhar ankowski; do
      "$GEN" --events "$NE" --channel "$ch" --model "$m" --seed "$SEED" --fsi off \
             2>> report/forbidden_log.txt > "data/win_${ch}_${m}.txt"
   done
done

echo ">> Rendering all figures with shared style..."
root -l -b -q macros/plot_proton.C
root -l -b -q macros/plot_kaon.C
root -l -b -q macros/plot_removal.C
root -l -b -q macros/plot_models.C
root -l -b -q macros/plot_kaon_models.C
root -l -b -q macros/compare_erem.C
root -l -b -q macros/compare_kaon_binding.C
root -l -b -q macros/compare_fsi.C
root -l -b -q macros/fsi_outcomes.C
root -l -b -q macros/fsi_spectra_channels.C
root -l -b -q macros/fsi_yield.C
root -l -b -q macros/fsi_multiplicity.C
root -l -b -q macros/window_eff.C
root -l -b -q macros/plot_spectral.C
root -l -b -q macros/plot_ankowski_sf.C
root -l -b -q macros/plot_xsec.C
root -l -b -q macros/offshell_W.C

# DUNE event-rate predictions: fold the MC efficiency tables (window + FSI) into
# the per-channel expected counts, then plot them. Must run after the macros that
# write report/{window_summary,fsi_summary}.txt above.
echo ">> Predicting DUNE event rates across all channels..."
./build/EventPredictor          # writes report/event_predictions.txt
root -l -b -q macros/plot_predictions.C

echo ">> Done. All figures in plots/ rebuilt with consistent styling:"
ls -1 plots/*.png
