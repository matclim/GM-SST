# SHiP straw tracker — simulation and reconstruction

A standalone Geant4 + GeoModel simulation of the SHiP spectrometer straw tracker,
with an ACTS-based reconstruction chain that fits tracks and reconstructs decay
vertices. It is built for two use cases: characterising the tracker with single
particles, and studying displaced multi-body decays (e.g. a long-lived particle
decaying to four pions).

Two executables:

- **`run_StrawTracker`** — the Geant4 simulation. Produces a hits file and,
  optionally, the straw geometry table.
- **`run_reco`** — the reconstruction. Reads the hits, fits tracks, and
  reconstructs vertices.

For the detailed physics and data-flow, see **[DOCUMENTATION.md](DOCUMENTATION.md)**.
This file is the quick start.

---

## Building

Requirements: Geant4 (≥ 11.4), GeoModel, ROOT, and ACTS (core, built from
`main`; this code targets the API around commit `b660c71`). Point CMake at your
ACTS install and build:

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build . -j
```

`RelWithDebInfo` is recommended: a plain `Debug` build re-enables an ACTS
field-interpolation assertion that fires on benign edge cases, and a pure
`Release` build has caused optimiser-exposed crashes in the past. `RelWithDebInfo`
keeps the assertions off and the symbols on.

---

## Coordinates

Two frames appear throughout, and mixing them up is the single most common source
of confusion:

- **SHiP frame** — the experiment's frame, in mm. The generator file, the output
  trees, and all physics discussion use this.
- **World frame** — the internal Geant4/GeoModel frame, `= SHiP − 60000 mm`. The
  straw geometry table and the `--field-origin-z` option are in this frame.

`run_reco` writes its output positions in the **SHiP frame** and prints the
station z in both, so you can always cross-check.

Key positions (SHiP frame, mm): the four stations at 84320 / 86820 / 91820 /
93320; the magnet centre at 89220; the decay volume spanning 33120 → 83120.

---

## Quick start

### 1. Dump the straw geometry (once per geometry)

```bash
./run_StrawTracker --dump-straws straws.root
```

This writes a `Straws` tree with every wire's position and direction, in the
**world** frame. It is the contract between sim and reco: the reco builds one
ACTS surface per row, so this file must come from the same geometry that produced
the hits.

### 2a. Single particles (tracker characterisation)

```bash
MAP=path/to/MainSpectrometerField.root

./run_StrawTracker --n-events 500 --particle mu- --energy-MeV 10000 \
  --pos-mm 0 0 79320 --dir 0 0 1 --field-map $MAP --output mu.root

./run_reco --hits mu_t-1.root --field $MAP --straws straws.root \
  --field-origin-z 29220 --kbl 0.84 --meas-res 0.16 \
  --select primary --p-seeds 2,5,10,15,20,30 --reco-out mu_reco.root
```

### 2b. K⁰_S → π⁺π⁻ (two-body vertexing)

```bash
./run_StrawTracker --n-events 500 --particle kaon0S --energy-MeV 10000 \
  --pos-mm 0 0 79320 --dir 0 0 1 --force-decay 2cpi \
  --field-map $MAP --output ks.root

./run_reco --hits ks_t-1.root --field $MAP --straws straws.root \
  --field-origin-z 29220 --kbl 0.84 --meas-res 0.16 \
  --select signal --p-seeds 1,2,3,5,7,10,15,20 --reco-out ks_reco.root
```

### 2c. Long-lived particle → 4π (displaced multi-body vertexing)

```bash
./run_StrawTracker --n-events 500 --llp-file path/to/DP_4pi.root \
  --field-map $MAP --output dp4pi.root

./run_reco --hits dp4pi_t-1.root --field $MAP --straws straws.root \
  --field-origin-z 29220 --kbl 0.84 --meas-res 0.16 --select llp \
  --p-seeds 1,2,3,5,7,10,15,20,30,50,80,120 --chi2-max 50 \
  --reco-out dpvtx.root
```

### 3. Look at the results

```bash
root -l dpvtx.root
root [0] Vertices->Draw("rz>>h(200,-2000,2000)","fitOK==1"); h->Fit("gaus","","",-500,500)
root [1] Tracks->Draw("chi2/nMeas")     # should peak near 1
```

Always fit the **core** with a Gaussian on a restricted range; the distributions
have heavy tails and the bare RMS is outlier-dominated.

---

## Output at a glance

`run_reco` writes three trees to the `--reco-out` file:

- **`Vertices`** — one row per event with ≥ 2 fitted tracks (successes and
  failures both, distinguished by `fitOK`): fitted position, truth, residual,
  seed, quality (DOCA, impact parameters), the reconstructed parent's impact
  parameter to the target (`ipToOrigin`) and its momentum, and acceptance counts.
- **`Tracks`** — one row per fitted track: fitted vs truth momentum, angles,
  charge, χ², and acceptance flags.
- **`Hits`** — per-straw diagnostics: the drift measurement vs truth, and the
  left/right assignment.

Full branch listings and physics notes are in **[DOCUMENTATION.md](DOCUMENTATION.md)**.

---

## Current status

The chain runs end to end for single particles, K⁰_S, and LLP → 4π. For the LLP
sample the vertex resolution is a clean funnel: σ_z ≈ 150 mm for decays near the
tracker, widening to ~800 mm for the most upstream decays, unbiased throughout.
The near-tracker resolution is close to the geometric limit set by the daughter
opening angle; the far-field width is dominated by the long extrapolation and is
largely irreducible for this geometry. See the documentation for the open items.
