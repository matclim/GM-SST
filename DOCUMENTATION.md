# SHiP straw tracker — reconstruction documentation

This document explains the reconstruction chain in enough detail that someone new
to the code can run it, read its output, and understand why each step is there.
It assumes familiarity with tracking concepts (helix fits, Kalman filtering) but
not with this codebase or with ACTS.

The companion [README.md](README.md) is the quick start; this is the reference.

---

## 1. Overview and data flow

```
              Geant4 + GeoModel                     ACTS
   ┌──────────────────────────────┐   ┌──────────────────────────────────┐
   │  run_StrawTracker             │   │  run_reco                         │
   │                               │   │                                   │
   │  primary generator            │   │  read hits + straw table          │
   │    · particle gun             │   │  build straw surfaces + field     │
   │    · LLP decay replay         │   │  per track:                       │
   │  Geant4 tracking through      │──▶│    seed → Kalman fit              │
   │    the straws (with material) │   │  per event:                       │
   │  straw hits + truth           │   │    perigee-reference tracks       │
   │                               │   │    Billoir vertex fit             │
   │  → hits file (ROOT)           │   │  → Vertices / Tracks / Hits trees │
   │  → straw geometry table       │   └──────────────────────────────────┘
   └──────────────────────────────┘
```

The two halves communicate through two files: the **hits file** (per-straw
energy deposits, truth information, and drift times) and the **straw geometry
table** (every wire's position and direction). The straw table is the contract:
the reco builds one ACTS surface per row, so it must be produced by the same
geometry that produced the hits.

---

## 2. Coordinates and geometry

Everything is quoted in the **SHiP frame** (mm) at the interfaces. Internally the
simulation works in a **world frame** offset by

```
world = SHiP − 60000 mm
```

| element | SHiP z (mm) | world z (mm) |
|---|---|---|
| decay volume | 33120 … 83120 | −26880 … 23120 |
| gun (default) | 79320 | 19320 |
| station 0 | 84320 | 24320 |
| station 1 | 86820 | 26820 |
| magnet / field-map centre | 89220 | 29220 |
| station 2 | 91820 | 31820 |
| station 3 | 93320 | 33320 |

Volume hierarchy (each contains the next):

```
World          ±6000 × ±6000 × ±60000 mm      SHiP z 0 … 120 m
└─ Spectrometer ±5000 × ±5000 × ±10000 mm     20 m, centred on the magnet
   │            x/y match the field map exactly; THE FIELD LIVES HERE
   ├─ Station_0 … Station_3                    aperture ±2000 (x) × ±3000 (y)
```

Each station: 4 layers × 2 sub-layers × 300 straws = 2400 straws (9600 total).
Straws are stereo (±2.3°, alternating by layer); sub-layers are staggered by one
radius. Straw radius 10 mm, length 4000 mm, pitch 20 mm.

**Why the Spectrometer volume exists.** The field map extends ±7500 mm in z, so
it overlaps the whole tracker. A separate magnet box would either overlap the
stations (a Geant4 error) or, if made smaller, silently truncate the field
integral — which under-bends the tracks and makes the reconstructed momentum too
high. Sizing one volume to the map and putting the stations inside it means the
field acts everywhere the map defines it, and a station moved outside the mapped
region fails loudly rather than reconstructing at the wrong momentum.

---

## 3. The magnetic field

The field comes from a ROOT map (`Range` tree: grid bounds and cell size; `Data`
tree: the field samples). **Its coordinates are in centimetres and its field
values in Tesla** — both `ShipFieldMap` (sim) and `ActsShipField` (reco) convert
the coordinates cm → mm on read. Reading them as mm shrinks the field region
tenfold, which un-bends any track more than half a metre off axis; the code
applies the conversion in one place each, and the loader prints its true extent
(±5000 mm in x/y, ±7500 mm in z) so a regression is visible.

The integral along the axis is ≈ 0.84 T·m, giving ~25 mrad of bend for a 10 GeV
track. `run_reco --dump-field` prints the field the reco actually sees, along its
z axis and transversely, plus the integral — use it whenever the field is in
doubt. A straight 10 GeV muon's deflection at station 3 is the end-to-end check:
a silently wrong field is invisible in track *directions* and *vertices* (which
come from the straws) and shows up only in the momentum.

---

## 4. Simulation — `run_StrawTracker`

Geant4 11.4 with GeoModel geometry and the FTFP_BERT physics list (includes decay
and hadronic physics, so pions can interact — about 6% are lost to nuclear
interactions across the four stations, which is real and appears as the gap
between geometric and hit acceptance in the reco output).

### Primary generators

**Particle gun** (default) — one particle, configured from the CLI:

```
--particle <name>     e.g. mu-, pi-, kaon0S
--energy-MeV <E>
--pos-mm <x> <y> <z>  SHiP frame
--dir <dx> <dy> <dz>
--force-decay 2cpi    force K0_S → π+ π− (installed after Initialize)
```

**LLP replay** (`--llp-file <file>`) — reads an external decay sample and fires
every **charged** daughter as a primary from the recorded decay vertex. The LLP
itself is never tracked (it is neutral and has already decayed).

**Neutral particles.** A straw tracker measures ionisation, so it is blind to
neutral particles by construction — a neutron, photon, or pi0 does not ionise the
gas and leaves no hits. The chain reflects this at two stages: the generator does
not even fire neutral daughters (it skips any daughter with zero PDG charge, since
they would leave no signal), and the reconstruction's `llp` selection keeps only a
charged allow-list (pi/K/p/mu/e). So neutrals are neither simulated nor
reconstructed — there is nothing to reconstruct. This is detector physics, not a
limitation of the code. Two consequences worth keeping in mind:

- *Invariant mass undershoots for modes with neutrals.* Only the charged daughters
  enter the summed four-momentum, so if the true final state includes a neutral
  (e.g. pi+ pi- pi0), the reconstructed `invMass` sits below the true parent mass,
  with a low-side tail. The pure DP -> 4pi sample is all-charged, so this does not
  arise there — but any mode with neutrals carries this bias, and it should be
  modelled. The vertex and pointing are unaffected: they are built from the
  charged tracks, which still converge at the true decay point.
- *Neutrals that interact* (a neutron scattering, a photon converting) produce
  charged secondaries, but those have `parentID != 0` and so are correctly
  excluded by the `parentID == 0` cut in the `llp` selection — they are not
  mistaken for LLP daughters. The input tree
carries the decay vertex (`vtx_x/y/z`, in metres, SHiP frame) and the daughters
(`d_px/py/pz/E`, in GeV, as vectors), converted m → mm and GeV → MeV on read.

### Other options

```
--field-map <file>    the ROOT field map
--output <file>       hits file (written as <stem>_t-1.root)
--dump-straws <file>  write the straw geometry table and exit
--n-events <N>
--seed <s>
```

### The hits file

A per-straw hit tree carries: the track and particle IDs (`trackID`,
`parentID`, `pdg`); the straw address (`stationID`, `layerID`, `subLayerID`,
`strawID`); the step entry/exit points; the **production vertex** (`vtxX/Y/Z`)
and **truth momentum at production** (`vpx/vpy/vpz`, MeV) — these make every
truth comparison downstream possible; and the **drift time** (`driftTime`, ns)
plus the true distance of closest approach (`driftTrue`, mm, diagnostic only).

### The drift model

A real straw does not measure a position — it measures the arrival **time** of
the first ionisation electrons at the wire. The simulation reproduces this:

1. a charged track leaves discrete primary-ionisation clusters along its path
   (Ar/CO₂: ~28 clusters/cm);
2. each cluster's electrons drift to the wire; the TDC records the earliest
   arrival, so the measured radius is set by the *nearest* cluster, not the
   track's true closest approach;
3. the arrival time is smeared by diffusion (growing with drift distance) and
   TDC jitter.

This makes the resolution **radius-dependent** — worst near the wire (few
clusters lie near the closest point, and the track runs tangential to the drift
circles) and near the wall (diffusion), best at mid-radius. The delivered
resolution is ~100–160 µm, following NA62's straws (same gas, same technology;
our straws are larger, so slightly worse). The model lives in `StrawDrift.h` and
is shared by sim and reco, so there is a single r-t relation.

The reco inverts the drift time to a radius with a **calibrated** r-t relation
(a polynomial fit to ⟨true DOCA | t⟩), because the nearest-cluster effect biases
the naive linear inversion outward near the wire — exactly as a real experiment
calibrates its r-t from data.

---

## 5. Reconstruction — `run_reco`

### 5.1 Inputs and options

```
--hits <file>            the hits file
--straws <file>          the straw geometry table
--field <file>           the field map
--field-origin-z <mm>    field-map centre in the WORLD frame (29220)
--kbl <T·m>              field integral for the seeder (0.84)
--meas-res <mm>          declared measurement resolution (0.16 matches the drift model)
--select <mode>          which tracks to reconstruct (see below)
--p-seeds <list>         momentum seed grid, GeV (comma-separated)
--chi2-max <x>           reject a track if its best χ²/ndf exceeds this (50)
--reco-out <file>        output trees
--n <N>                  process only the first N events
--ship-z-origin <mm>     output-frame offset (60000; positions written in SHiP frame)
--ip-origin-z <mm>       target z (SHiP) for the reconstructed-parent impact parameter (0)
--daughter-mass <GeV>    mass hypothesis for the invariant mass (0.13957 = pion)
--truth-pid              use each daughter's TRUE pdg for its mass (perfect-PID reference)
--n-events <N>           process only N events (alias of --n; default all)
--start-at-event <K>     skip to event index K before processing (default 0)
--dump-field             print the reco's field and exit
--verbose-fail <N>       print the first N propagation failures in detail
```

### 5.2 Track selection (`--select`)

| mode | keeps | use |
|---|---|---|
| `signal` | `parentID==1 && \|pdg\|==211` | K⁰_S → π⁺π⁻ daughters |
| `llp` | `parentID==0` and charged | LLP daughters (fired as primaries) |
| `primary` | `trackID==1` | single-particle gun |
| `all` | everything | debugging |

### 5.3 Geometry and field

`ShipStrawGeometry` reads the straw table and builds one ACTS `StrawSurface` per
wire, keyed by a `GeometryIdentifier` encoding (station, layer, sub-layer,
straw). Each surface carries a homogeneous material slab (an effective Mylar
thickness matching the real straw's radiation-length budget), so the Kalman
filter accounts for multiple scattering. `ActsShipField` builds the interpolated
field from the same ROOT map.

### 5.4 The measurement, and left/right

For each straw a track crossed, the reco takes the earliest drift time across the
Geant4 steps in that straw and inverts it to an **unsigned** radius. A straw
gives no information about *which side* of the wire the track passed — the
left/right ambiguity — so the sign is not supplied by truth; it is resolved by
the fit. The calibrator takes the side from the Kalman filter's predicted state
at each surface.

### 5.5 Seeding

The sagitta-based momentum estimate is too noisy to trust, so momentum and charge
are found by a **scan**: for each track, and each (momentum seed, charge sign) on
the grid, the track is fitted and the best χ²/ndf kept; the track is rejected if
even the best exceeds `--chi2-max`. Two subtleties matter and are handled:

- the q/p covariance is rebuilt from the grid momentum, not reused from the
  seeder — otherwise the fit is pinned to the seed and cannot move;
- a track is *rejected*, not clamped, if no seed fits — a fit allowed to run to
  very high momentum finds a nearly-straight track that fits either charge
  equally well, which destroys the charge determination.

### 5.6 Track fit

An ACTS `KalmanFitter` with a `DirectNavigator` over the per-event sequence of
straw surfaces (sorted in z), multiple scattering and energy loss enabled. The
reference surface is a plane just upstream of the first straw.

### 5.7 Vertex fit

Works for any N ≥ 2 tracks:

1. **Seed** — the least-squares point closest to all track lines.
2. **Perigee-reference** — each track is re-expressed at its point of closest
   approach to the seed vertex, using ACTS's `ImpactPointEstimator`. This is
   essential: Billoir expects perigee parameters *with respect to the vertex* and
   linearises around it. For a displaced decay this can be a back-extrapolation
   of tens of metres, which a naive "propagate to a line at the seed" cannot do
   reliably — the estimator finds the closest-approach point properly.
3. **Fit** — a `FullBilloirVertexFitter` with a `HelicalTrackLinearizer`. The
   seed sets the linearisation point only, so the result is a genuine fit, not a
   pull toward the seed.
4. **Quality** — pairwise track-to-track DOCA and each track's impact parameter
   with respect to the fitted vertex.

---

## 6. Output trees

Positions are written in the **SHiP frame**; residuals and covariances are
frame-independent.

### `Vertices` — one row per event with ≥ 2 fitted tracks

Rows are written for **failures** too (see `fitOK` / `vtxFail`), so inefficiency
is visible rather than silently absent.

| group | branches | meaning |
|---|---|---|
| identity | `event`, `nTrk` | event id, tracks in the fit |
| fitted | `vx vy vz` | fitted vertex (SHiP mm) |
| covariance | `sx sy sz` | position uncertainties |
| truth | `tx ty tz` | true decay vertex (SHiP mm) |
| residual | `rx ry rz` | fitted − truth |
| momentum/charge | `p0 p1 q0 q1` | leading two tracks |
| quality | `docaMax docaMean ipMax ipMean` | vertex consistency (see note) |
| pointing | `ipToOrigin ipCApZ` | reconstructed-parent impact parameter to the target |
| parent | `parentPx parentPy parentPz` | reconstructed parent momentum (GeV) |
| mass | `invMass` | invariant mass of the daughters (GeV, see note) |
| weight | `weight` | event weight (from the LLP file; 1 otherwise) |
| seed | `sdx sdy sdz seedRz seedCond` | seed position, its z residual, conditioning |
| status | `fitOK vtxFail nProp` | success flag; failure mode; tracks reaching the perigee |
| acceptance | `nTruthAcc nTruthHit nFitted`, `nGeoAcc0-3`, `nHitAcc0-3` | truth-in-acceptance, hit, fitted counts; per-station |

`vtxFail`: 0 = ok, 1 = fewer than 2 tracks, 2 = perigee propagation failed,
3 = Billoir failed.

**Two different "impact parameters" live in this tree — do not confuse them.**

- **`ipMax` / `ipMean`** measure *vertex quality*. For each track in the vertex,
  its impact parameter is taken *with respect to the fitted vertex itself* — how
  far that track passes from the reconstructed decay point. `ipMean` is the
  average over the tracks (overall consistency); `ipMax` is the worst single
  track (an outlier catcher: a mis-fit, wrong-side, or mis-associated track shows
  up here). Small values mean the daughters genuinely converge. Cut on `ipMax`
  (or `docaMax`) to reject poorly reconstructed vertices.

- **`ipToOrigin`** measures *pointing*, to a different reference. It is the impact
  parameter of the reconstructed **parent** — the line through the fitted vertex
  along the summed daughter momentum — to the **target** at (0, 0, `--ip-origin-z`)
  in the SHiP frame (default the SHiP origin). It answers a physics question, not
  a quality one: does the reconstructed decay point back to where the parent was
  produced? A prompt background tends to point back (small `ipToOrigin`); a
  genuine displaced decay has a characteristic distribution broadened by the
  direction resolution. `ipCApZ` is the SHiP z of the closest-approach point along
  the parent line — it tells you *where* the parent points, not just how far it
  misses. `ipToOrigin = −1` flags an undefined value (failed fit or zero parent
  momentum); cut `ipToOrigin >= 0` in plots.

  `parentPx/Py/Pz` are the reconstructed parent momentum (the vector sum of the
  fitted daughters, GeV) — useful in its own right, and the direction that
  defines the pointing line above.

### `Tracks` — one row per fitted track

| group | branches |
|---|---|
| identity | `event`, `pdg` |
| momentum | `pFit pTrue pRes` |
| direction | `phiFit phiTrue dphi thetaFit thetaTrue dtheta dangle` |
| charge | `qFit qTrue` |
| fit quality | `chi2 nMeas bestChi2` |
| acceptance | `inAcc`, `geoAcc0-3`, `hitAcc0-3`, `nGeoAcc`, `nHitAcc` |

`geoAcc` is whether a straight-line truth track crosses a station's aperture;
`hitAcc` is whether it actually left hits. The difference is the material loss.

### `Hits` — per-straw diagnostics (truth-level, for debugging)

`event trk station layer sub straw z`, then the measurement diagnostics:
`driftMeas` (the unsigned radius the fit used), `truthLoc0` (signed truth drift),
`wrongSide` (1 if the fit chose the wrong side), `driftTrue`, `diff`, `theta`,
`pTrue`, `chordDot`, `pathLen`, `nSteps`.

---

## 7. Analysis notes

**Always fit the core, never quote the bare RMS.** The distributions have heavy
tails (wrong-side hits, scattered or interacting pions, far-field vertices), and
the RMS is outlier-dominated. Fit a Gaussian on a restricted range:

```cpp
Vertices->Draw("rz>>h(200,-2000,2000)","fitOK==1");  h->Fit("gaus","","",-500,500);
Tracks->Draw("dtheta>>hc(200,-0.01,0.01)");          hc->Fit("gaus","","",-0.002,0.002);
Tracks->Draw("chi2/nMeas");                          // should peak near 1
Tracks->Draw("qFit==qTrue");                         // charge accuracy
```

**The vertex-resolution funnel** is the headline result for a displaced sample.
σ_z grows with the decay's distance from the tracker, because the vertex is found
by extrapolating the tracks back:

```
σ_z ≈ σ_θ · L / θ_open
```

with σ_θ the per-track angular resolution, L the extrapolation distance, and
θ_open the daughter opening angle. Plot it directly:

```cpp
Vertices->Draw("rz:tz","fitOK==1 && abs(rz)<2000","COLZ");
```

For the LLP → 4π sample: σ_z ≈ 150 mm near the tracker, ~800 mm for the deepest
decays. The near-tracker value is close to the geometric limit set by θ_open and
σ_θ; the far-field width is dominated by L and is largely irreducible for this
geometry.

**Pointing.** `ipToOrigin` (see the `Vertices` note in §6) is the reconstructed
parent's impact parameter to the target — a displaced-decay discriminant. Plot it
against decay z to see how the direction resolution feeds into pointing:

```cpp
Vertices->Draw("ipToOrigin:tz","fitOK==1 && ipToOrigin>=0","COLZ");
```

Expect the far decays, whose parent direction is least well measured, to have the
larger impact parameters — the same σ_θ that widens the vertex funnel also blurs
the pointing.

**χ²/ndf ≈ 1** indicates the measurement error model is honest — the declared
`--meas-res` matches the drift model's real spread (~160 µm). A value well below 1
means the fit is under-weighting the hits (declared error too large); well above
1 means either the error is too small or hits are being mis-assigned (e.g. wrong
left/right side).

---

## 8. Validation checks worth keeping

- **Field integral** — a straight 10 GeV muon's y-deflection at station 3. A
  silently misplaced or truncated field is otherwise invisible; only the momentum
  reveals it.
- **Straw table** — `Straws->Draw("cz")` must show the four station z's in the
  world frame; `run_reco` prints the z it derives, in both frames, every run.
- **χ²/ndf near 1** on a clean single-particle sample, confirming the measurement
  model is weighted correctly.

---

## 9. Open items and known limitations

- **Far-field vertex resolution** is dominated by the long back-extrapolation of
  near-parallel, stiff (50–120 GeV) daughters, whose longitudinal crossing is
  weakly constrained. This is close to the geometric limit; seeding tricks (a
  z-scan minimising transverse spread) were tried and did not improve on the
  least-squares seed, consistent with the constraint being intrinsic rather than
  algorithmic.
- **Left/right assignment** is resolved from the Kalman prediction. The
  `wrongSide` branch lets you measure the mis-assignment rate; a full Combinatorial
  Kalman Filter would instead branch on both hypotheses and let χ² decide.
- **Drift-model parameters** (cluster density, drift velocity, diffusion) follow
  NA62; if SHiP publishes measured values, re-fit the calibrated r-t in
  `StrawDrift.h` and update `--meas-res` to the resulting resolution.
- **Track finding** is not yet combinatorial — the reco groups hits by truth
  track ID. A CKF would be the next step for genuine pattern recognition.
- **Neutral final-state particles** are invisible to the tracker by construction
  (no ionisation, no hits) and are neither simulated nor reconstructed. For decay
  modes containing neutrals, the reconstructed invariant mass is biased low by the
  missing neutral momentum; this must be modelled per mode. All-charged modes
  (like DP -> 4pi) are unaffected.
