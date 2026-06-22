# Tabulated spectral functions S(p, E)

Grid spectral functions used by the `benhar` momentum model (`--model benhar`),
loaded by `include/PDKSpectral.h`. (The `ankowski` model builds its effective
S(p,E) analytically and does not read these grids.)

| file | nucleon |
|------|---------|
| `gsf_Ar40P.grid` | argon-40 proton |
| `gsf_Ar40N.grid` | argon-40 neutron |

**Source:** the NuWro generator, `data/sf/` of
<https://github.com/NuWro/nuwro> (raw e.g.
`https://raw.githubusercontent.com/NuWro/nuwro/master/data/sf/gsf_Ar40P.grid`).
These are the JLab E12-14-012 argon spectral functions
(arXiv:2203.01748, arXiv:2312.13369).

**Format** (NuWro `gridfun2d`), units **MeV**:

```
eRes pRes          # grid counts: #energy points, #momentum points (200 200)
eMin pMin          # axis minima (0 0)
eMax pMax          # axis maxima (400 800)  -> E in [0,400] MeV, p in [0,800] MeV
# then pRes momentum blocks, each:
p                  # momentum [MeV]
e  S(p,e)          # eRes (energy[MeV], value) pairs (wrapped, several per line)
...
```

Values are stored momentum-outer / energy-inner. Only relative `S` values matter
for sampling; the loader converts the axes from MeV to GeV.
