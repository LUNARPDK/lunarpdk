# DUNE PDK Monte Carlo

Toy Monte Carlo for **bound-nucleon decay** (proton and neutron) in a DUNE 10 kt
liquid-argon module. The default channel is **p → K⁺ν̄**; a range of SUSY/GUT
two-body modes is supported (see [Decay channels](#decay-channels)).

## Layout

```
include/   Shared headers
  PDKChannels.h    Nucleon type, particle masses, decay-channel registry
  PDKConfig.h      Config struct, params.dat loader, polynomial sampler, RNG
  PDKKinematics.h  Off-shell two-body decay N -> d1 + d2 + Lorentz boost
  PDKMomentum.h    Nucleon momentum models (Fermi gas, LFG, SRC, spectral fn.)
  PDKSpectral.h    Tabulated S(p,E) loader + 2D sampler (benhar)
  PDKFsiXsec.h     Hadron-nucleon cross sections (Metropolis piN/NN + K/eta)
  PDKCascade.h     Intranuclear FSI cascade for the final-state hadrons
src/       Programs
  EventPredictor.cpp        Expected event count vs. nucleon lifetime
  PDKMCGenerator.cpp        Sanity check of the momentum sampler
  LunarPDKGenerator.cpp  Full event generator (writes the ascii table)
macros/    ROOT macros: load the ascii output into a TTree/TNtuple, and plot
  pdk_style.h      Shared plot style + model/binding colour palette (all macros)
  plot_spectral.C  S(p,E) heatmap of the proton/neutron NuWro grids (benhar)
  plot_ankowski_sf.C  S(p,E) heatmap of the analytic Ankowski effective SF
  plot_xsec.C      Hadron-nucleon FSI cross sections (includes PDKFsiXsec.h)
  fsi_multiplicity.C  Mean post-FSI final-state multiplicity per channel
  offshell_W.C     Off-shell mass W + kinematically-forbidden fraction by model
config/    params.dat — Fermi-momentum PDF (polynomial coeffs + range + f_max)
  fsi/             NuWro Metropolis NN cross-section tables (provenance only)
data/      Generator output (output.txt) and ROOT files (output.root)
plots/     Generated figures (proton_p.png, kaon_p.png, proton_models.png,
           spectral_pe.png, ankowski_sf.png, fsi_xsec.png, fsi_multiplicity.png,
           offshell_W.png)
build/     Compiled binaries
run.sh            One-command pipeline driver (single model)
compare_models.sh Overlay nucleon momentum across all nuclear models
analyze_kaon.sh   Daughter spectrum across all models (+ summary table)
compare_fsi.sh    Overlay the meson spectrum with vs. without FSI (kaon + pion)
compare_fz.sh     Primary-meson outcome with the formation zone off vs. on
make_plots.sh     Rebuild every figure in plots/ with the shared styling
```

All ROOT macros include `macros/pdk_style.h`, which sets one global style and a
fixed colour/label map so a given model (or binding) keeps the same colour and
name in every figure. To regenerate the whole figure set in one go:

```sh
./make_plots.sh [n_events] [seed] [channel]   # defaults: 200000 events, seed 1, pToKnu
```

## Quick start (from the project root)

```sh
./run.sh [n_events] [seed] [model] [binding] [channel]   # defaults: 10000 1 polynomial potential pToKnu
```

This builds the binaries, generates events, fills the TTree
(`data/output.root`), and writes the momentum plots `plots/proton_p.png`
(initial nucleon momentum) and `plots/kaon_p.png` (hadron-side daughter).

To compare the proton momentum distributions of all nuclear models on one plot:

```sh
./compare_models.sh [n_events] [seed] [binding]   # writes plots/proton_models.png
```

To analyze the lab-frame **kaon** momentum spectrum (the observable p → K⁺ν
signature) across all models:

```sh
./analyze_kaon.sh [n_events] [seed] [binding]
# writes plots/kaon_models.png and report/kaon_summary.txt
```

A free proton at rest emits a monochromatic kaon at p_K = (M_p²−m_K²)/(2M_p) ≈
0.339 GeV/c; Fermi motion broadens it and the binding (off-shell W < M_p) shifts
the peak down. The summary table reports, per model, the mean/RMS/peak kaon
momentum, the fraction inside a nominal [0.30, 0.38] GeV/c signal window, and the
low-momentum tail fraction (p_K < 0.25 GeV/c) fed by deep binding and SRC pairs.

To isolate the **binding-model** dependence of the kaon spectrum (potential vs.
constant vs. shell), for two analytic momentum models (`lfg`, `sf`) — binding
does not apply to the spectral-function models `benhar`/`ankowski`:

```sh
for b in potential constant shell; do for m in lfg sf; do
  ./build/LunarPDKGenerator --events 300000 --model "$m" --binding "$b" --seed 7 \
    > "data/kbind_${m}_${b}.txt"
done; done
root -l -b -q macros/compare_kaon_binding.C
# writes plots/kaon_binding.png and report/kaon_binding_summary.txt
```

The binding model sets the removal energy, hence the off-shell mass W < M_p,
hence the kaon momentum. The result is that the binding choice mostly produces a
rigid ~5–10 MeV/c shift of the peak position — the momentum-dependent optical
potential is the deepest and gives the softest kaon — while the *width* of the
spectrum is governed by the momentum model, not the binding (see the report).

## Nucleon momentum models

Selected with the `--model` flag. All are toy/parametrized approximations for
argon-40 (k_F ≈ 0.217 GeV/c protons, ≈ 0.230 GeV/c neutrons), defined in
`include/PDKMomentum.h`.

| model        | description |
|--------------|-------------|
| `polynomial` | toy PDF read from `config/params.dat` (default) |
| `gfg`        | global Fermi gas: n(p) ∝ p² for p < k_F |
| `lfg`        | local Fermi gas: Fermi spheres with k_F(r) from a Woods-Saxon density |
| `src`        | Fermi-gas bulk + a 1/p⁴ short-range-correlation tail above k_F |
| `sf`         | toy spectral function: local-Fermi-gas bulk + SRC tail (parametrized) |
| `benhar`     | **tabulated** S(p,E) from a NuWro grid file (default: real JLab argon SF) — draws (p, E_rem) jointly from the data |
| `ankowski`   | **analytic** effective Ankowski-Sobczyk S(p,E): shell mean field (HO momentum + shell separation energy) + a correlated 1/p⁴ tail |

The `benhar` model reads a real tabulated spectral function S(p, E)
(`include/PDKSpectral.h`) instead of a parametrization. The grid file is chosen by
parent nucleon — `config/sf/gsf_Ar40P.grid` (proton) or `gsf_Ar40N.grid`
(neutron) — or set explicitly with `--sf-file PATH`. See `config/sf/README.md` for
the data source (NuWro / JLab E12-14-012) and the grid format.

The `ankowski` model builds an effective spectral function analytically (no grid):
each shell contributes its harmonic-oscillator momentum profile paired with its
measured separation energy, and a correlated fraction (`as_corr_fraction`) carries
a 1/p⁴ tail with the two-nucleon removal energy E ≈ E_offset + p²/2M.

Because both spectral-function models already encode the full momentum **and**
removal-energy distribution, `--binding` does **not** apply to either of them.

## Binding-energy (removal-energy) options

The removal (separation) energy of a *mean-field* nucleon is set by a separate,
independently selectable binding model. Correlated SRC-tail nucleons always carry
the two-nucleon removal energy E ≈ E_offset + p²/2M regardless of this choice.

| binding     | description |
|-------------|-------------|
| `potential` | momentum-dependent optical potential V(k_F, p) [nucl-th/0311051] (default) |
| `constant`  | fixed average separation energy (`e_sep_const`, 30 MeV by default) |
| `shell`     | argon-40 nucleon shell levels, Gaussian-smeared about their separation energies (`argon_proton_shells()` / `argon_neutron_shells()`) |

The binding model applies to every analytic momentum model (all except `benhar`
and `ankowski`, which carry their own removal energy), e.g.
`./run.sh 10000 1 lfg shell` or
`./build/LunarPDKGenerator --model gfg --binding constant`.

## Decay channels

Select the mode with `--channel` (default `pToKnu`). All are two-body
parent → (lepton-side) + (hadron-side) decays; the lepton-side ν̄ is invisible.
Run `./build/LunarPDKGenerator --help` to list them.

| key | mode | | key | mode |
|-----|------|-|-----|------|
| `pToKnu`   | p → K⁺ ν̄  | | `nToEPim`  | n → e⁺ π⁻ |
| `pToEPi0`  | p → e⁺ π⁰ | | `nToMuPim` | n → μ⁺ π⁻ |
| `pToMuPi0` | p → μ⁺ π⁰ | | `nToNuPi0` | n → ν̄ π⁰ |
| `pToNuPip` | p → ν̄ π⁺  | | `nToNuEta` | n → ν̄ η |
| `pToEEta`  | p → e⁺ η  | | `nToNuK0`  | n → ν̄ K⁰ |
| `pToMuEta` | p → μ⁺ η  | | `nToEKm`   | n → e⁺ K⁻ |
| `pToEK0`   | p → e⁺ K⁰ | | | |
| `pToMuK0`  | p → μ⁺ K⁰ | | | |

Particle/nucleon masses and the channel registry live in
`include/PDKChannels.h`. Neutron channels use the neutron mass, a neutron Fermi
fraction (N/A) and a neutron shell table (`argon_neutron_shells()`).

## Final-state interactions (FSI)

The hadron produced in the decay (K⁺, K⁰, π, η, …) does not escape the argon
nucleus untouched: it re-interacts on the way out. A semi-classical
**intranuclear cascade** transports it (and every secondary it makes) through
the residual nucleus before the final state is reported. It is **on by default**;
turn it off with `--fsi off` (which restores the legacy momentum-only columns).

```sh
./build/LunarPDKGenerator --channel pToEPi0 --fsi on --events 100000   # default
./build/LunarPDKGenerator --channel pToKnu  --fsi off                  # legacy output
```

| flag         | meaning |
|--------------|---------|
| `--fsi on\|off` | run (default) or skip the hadron FSI cascade |
| `--fsi-pot`     | apply a constant nucleon exit potential (~40 MeV) at the surface |
| `--formation-zone` | free-stream each produced meson one formation length before FSI is enabled |
| `--formation-time T` | formation time c·τ_f in fm (default 0.342; only with `--formation-zone`) |
| `--decay-mesons`| decay the escaped unstable mesons (π⁰/η → γγ, K⁰ → K_S/K_L, K_S → ππ); FSI-on only |

**Model** (`include/PDKCascade.h`, following NuWro's `kaskada`): the decay vertex
sits where the decaying nucleon is — the local-density models (`lfg`, and the
mean-field part of `sf`) already sampled a radius when they drew the nucleon, and
the cascade reuses it, sampling only the direction; models with no radius of
their own draw one independently from the same r²ρ(r) weight. The hadron is then
stepped in 0.05 fm
segments with a local mean free path λ = 1/(ρ·σ); on interaction a target
nucleon is drawn from a local Fermi sphere of radius k_F(r) — the spectators are
always a local Fermi gas, independently of the `--model` chosen for the decaying
nucleon — a sub-channel
(elastic / charge-exchange / absorption / production) is chosen, the kinematics
are generated in the projectile+target rest frame, and Pauli blocking is applied
to every outgoing nucleon. Secondaries are cascaded to completion. Four-momentum
and charge are conserved exactly in every channel. Leptons (e⁺/μ⁺/ν̄) do not
interact strongly and pass through.

**Cross sections** (`include/PDKFsiXsec.h`): the πN and NN cross sections are the
Metropolis tables (Phys. Rev. 110 (1958)) used by NuWro, embedded verbatim (the
raw NuWro files are kept under `config/fsi/` for provenance). The K⁺/K⁰ (S=+1,
nearly transparent), K̄ (S=−1, strongly absorbed via K̄N→Yπ) cross sections are
short interpolation tables of the measured K-nucleon totals (Dover & Walker 1982;
Friedman & Gal 2007; PDG), and the η (through the N\*(1535)) is an N\*(1535)
Breit–Wigner. A neutral kaon is propagated as the strangeness eigenstate it is
produced in: the nucleon-decay modes emit a K⁰ (S = +1), which is nearly
transparent like the K⁺, while a K̄⁰ (S = −1) made by charge exchange is strongly
absorbed. The strong interaction conserves strangeness, and the nuclear traversal
(~fm/c ~ 1e-23 s) is twelve orders of magnitude shorter than the K⁰_S lifetime,
so no K⁰–K̄⁰ mixing develops inside the nucleus; mixing and decay are left to the
out-of-nucleus stage in `PDKDecay.h`.
Quasi-elastic and charge-exchange scattering is forward-peaked (`dσ/dt ∝ exp(B·t)`),
with representative slopes B = 4, 3, 5 GeV⁻² for πN, KN and NN, reducing to
isotropic at threshold.

**Formation zone** (`--formation-zone`, off by default): a freshly produced
hadron is not yet a fully formed colour singlet and cannot re-interact at full
cross section immediately. With the flag on, each produced meson free-streams a
formation length `L_f = (p/m)·c·τ_f` before FSI is enabled, where `p/m = βγ`
dilates the rest-frame formation time `c·τ_f` (default 0.342 fm, SKAT;
`--formation-time` overrides) to the lab frame. The suppression applies to the
primary meson and to mesons made in inelastic production; struck recoil nucleons
are unaffected. The effect raises the escape (`none`) fraction — modestly for the
slow, nearly-transparent K⁺, strongly for the heavily-absorbed π⁰.

Each event is labelled by the **fate of the primary meson**: `none` (escaped
untouched), `elastic`, `cex` (left as a different meson), `produced` (survived but
made extra mesons), or `absorbed` (did not escape). Representative behaviour in
argon: the K⁺ is nearly transparent (~99% survives, almost all `none`/`elastic`,
<1% absorbed), while the π⁰ is heavily reworked (~39% `absorbed`, ~11% `cex`,
~13% `produced`); the K̄ and η are strongly absorbed.

This is a compact model for spectrum/efficiency studies: quasi-elastic and
charge-exchange scattering is forward-peaked, while absorption and multi-body
production are isotropic in the CM and use an approximate sequential phase space.

To visualize the effect, overlay the meson spectrum with vs. without FSI for the
transparent (K⁺) and strongly-interacting (π⁰) benchmark channels:

```sh
./compare_fsi.sh [n_events] [seed] [model]   # defaults: 200000 1 lfg
# writes plots/fsi_compare.png
```

Each curve is normalized per generated decay, so the FSI-on histogram shows both
the **depletion** (smaller area = mesons lost to absorption / charge exchange)
and the **softening** (a low-momentum tail from quasi-elastic energy loss). The
K⁺ stays ~99% transparent with a mild shift; the π⁰ drops to ~60% yield with a
pronounced low-momentum bump — the canonical "pion FSI eats the signal" effect.

## Running the steps manually

```sh
make
# The ROOT macros below read the legacy momentum-magnitude columns, so generate
# with --fsi off. For the post-FSI per-particle final state, drop --fsi off.
./build/LunarPDKGenerator --events 10000 --channel pToEPi0 --model lfg --seed 42 --fsi off > data/output.txt
root -l -b -q macros/text2tree.C    # ascii -> TTree
root -l -b -q macros/plot_proton.C  # TTree -> plots/proton_p.png (initial nucleon)
root -l -b -q macros/plot_kaon.C    # TTree -> plots/kaon_p.png   (hadron daughter)
```

`LunarPDKGenerator [--events N] [--channel KEY] [--model NAME] [--binding NAME]
[--fsi on|off] [--fsi-pot] [--seed S] [--config PATH]` — all flags are optional;
passing `--seed` makes the run reproducible, omitting it draws a non-deterministic
stream.

## Output columns (`data/output.txt`)

**With `--fsi off`** (legacy, momentum-magnitude table):

| column      | meaning                                  |
|-------------|------------------------------------------|
| `event`     | 1-based event index                      |
| `nucleon_p` | parent nucleon momentum magnitude [GeV/c] |
| `d1_p`      | lepton-side daughter momentum, lab frame [GeV/c] (e⁺/μ⁺/ν̄) |
| `d2_p`      | hadron-side daughter momentum, lab frame [GeV/c] (meson) |
| `e_rem`     | removal (separation) energy of the bound nucleon [GeV] |

**With `--fsi on`** (default, full post-FSI final state): one row per outgoing
particle, preceded by a per-event `# event <i>: nucleon_p=… e_rem=… outcome=…`
comment. Columns: `event  pdg  px  py  pz  E  outcome`, with momenta/energy in
GeV and `outcome` the primary-meson fate (`none`/`elastic`/`cex`/`produced`/
`absorbed`). Each event emits the visible lepton (PDG −11/−13) plus every hadron
that escaped the nucleus (PDG: π⁺ 211, π⁰ 111, π⁻ −211, K⁺ 321, K⁰ ±311,
K⁻ −321, η 221, p 2212, n 2112, Λ 3122). The particle count varies per event —
absorption removes the meson, knockout/production add nucleons and pions.

## Notes / approximations

- The bound nucleon is off-shell: mean-field nucleons are bound by the
  momentum-dependent optical potential V(k_F, p) of Juszczak, Nowak & Sobczyk,
  Eur. Phys. J. C 39 (2005) 195 [nucl-th/0311051], Eq. (8). Its nuclear-frame
  energy is E = sqrt(p^2 + M_N^2) + V(k_F, p), which sets the removal energy
  and the off-shell invariant mass W used for the decay. Deeply-bound,
  high-momentum nucleons with W < m₁ + m₂ are kinematically forbidden and
  resampled.
- The two-body decay N → d1 + d2 uses the general off-shell kinematics
  (`decay_two_body` in `include/PDKKinematics.h`); the daughter rest-frame
  momentum is p* = √λ(W², m₁², m₂²)/(2W). Unstable mesons (π⁰, η, K⁰) are
  reported as-is by default; pass `--decay-mesons` to decay the escaped ones
  (π⁰/η → γγ, K⁰ → K_S/K_L, K_S → ππ, recursing into the products).
- Final-state interactions of the hadron are modelled by the intranuclear
  cascade in `include/PDKCascade.h` (on by default, `--fsi off` to disable); see
  the [Final-state interactions](#final-state-interactions-fsi) section. It is a
  semi-classical spectrum-level model with Metropolis πN/NN cross sections,
  data-driven K-nucleon tables and an N\*(1535) Breit–Wigner for η,
  forward-peaked quasi-elastic/charge-exchange angular distributions and an
  approximate sequential phase space — not a precision transport code. Hyperons
  produced by K̄ absorption are emitted but not themselves cascaded or decayed.
- Neutron channels reuse the same momentum/binding machinery with the neutron
  mass, a neutron Fermi momentum/fraction (N/A) and a neutron shell table. These
  neutron nuclear inputs are representative argon values, not a fit; the small
  proton/neutron isospin asymmetry is otherwise approximate.
- The `gfg`/`lfg`/`src`/`sf` momentum models are parametrized toy approximations
  for argon, not tabulated nuclear inputs. The argon density (two-parameter
  Fermi, Eq. (4) of the same paper) and local Fermi momentum (Eq. (5)) feed both
  the local-density sampling and the binding potential. The `sf` model is a
  Fermi-gas + 1/p⁴ SRC-tail stand-in.
- The `benhar` model is **not** parametrized: it samples (p, E_rem) jointly from a
  real tabulated spectral function S(p, E) read from a NuWro grid file
  (`config/sf/gsf_Ar40{P,N}.grid`, the JLab E12-14-012 argon SF), overridable with
  `--sf-file`.
- The `ankowski` model is an **effective** spectral function built analytically
  (no grid): each argon shell contributes its harmonic-oscillator momentum profile
  paired with its measured separation energy, and a correlated fraction
  (`as_corr_fraction`) adds a 1/p⁴ tail with the two-nucleon removal energy. So its
  removal-energy spectrum shows the shell structure plus a correlated tail in the
  Ankowski-Sobczyk spirit, distinct from the `benhar` grid.
- The `shell` binding levels live in `argon_proton_shells()` /
  `argon_neutron_shells()` in `include/PDKMomentum.h`; their separation energies
  are representative Ar-40 values, not a fit.
- `EventPredictor` uses N = Np·(1 − e^(−T/τ)) ≈ Np·T/τ (valid for τ ≫ T).
