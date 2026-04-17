# StrawTracker — GeoModel + Geant4 Geometry

A full geometry description and Geant4 simulation of a straw tube tracker,
built with the GeoModel geometry toolkit.

---

## Detector layout

```
Beam direction: +Z

Station positions (z-centres):
  Station 0:  z = 26500 mm
  Station 1:  z = 29000 mm
  Station 2:  z = 34000 mm
  Station 3:  z = 35500 mm

Each station: 4 m (X) × 6 m (Y) active face
  └── 4 straw layers  (= 4 "views", stereo angle alternates ±2.3°)
        └── Material frame (one per view, hollow rectangle)
        └── 2 sub-layers per layer (staggered by ½ pitch)
              └── 300 straws per sub-layer (pitch = 2 cm)

Straw specification:
  Diameter  : 20 mm   (radius 10 mm)
  Length    : 4000 mm (along X, horizontal)
  Wall      : 30 µm Mylar (PET, ρ = 1.39 g/cm³)
  Fill      : Ar/CO₂ 70/30 by mass (ρ ≈ 1.56×10⁻³ g/cm³)

Stereo angles (rotation of layer volume around beam axis Z):
  Layer 0:  +2.3°   (u view)
  Layer 1:  −2.3°   (v view)
  Layer 2:  +2.3°   (u view)
  Layer 3:  −2.3°   (v view)

Staggering: sub-layer 1 is displaced by +½ pitch (+10 mm) in Y
            and +1 straw radius (+10 mm) in Z relative to sub-layer 0.

Frame per view (FairShip-style):
  Shape     : outer rectangle − inner aperture (GeoShapeSubtraction)
  Aperture  : 4020 mm (X) × 6030 mm (Y) — slight clearance beyond straw pattern
  Outer     : aperture + 100 mm bar on each side → 4220 × 6230 mm
  Thickness : 44 mm along Z (covers both sub-layers)
  Material  : configurable (default Aluminum)
```

---

## Volume hierarchy

```
World  (air box, ~15 m in Z)
 └── Station_i  [i = 0..3]   air box, placed at station z-centre
      └── Layer_j  [j = 0..3]   air box, rotated by ±2.3° around Z
           ├── StrawFrameLV_*   GeoShapeSubtraction  material frame (Aluminum)
           ├── SubLayer_0       air slab at z = −10 mm
           │    └── Straw_n     [n = 0..299]
           │         ├── StrawWall  GeoTube  rMin=9.97  rMax=10  Mylar
           │         └── StrawGas   GeoTube  rMin=0     rMax=9.97  ArCO₂  ← SD
           └── SubLayer_1       air slab at z = +10 mm, staggered
                └── Straw_n     [n = 0..299]
                     ├── StrawWall   …
                     └── StrawGas    …
```

---

## Magnetic field

A field volume (semi-transparent blue box, 4 m × 6 m × 3.5 m) is placed
between stations 1 and 2 with its lab-frame centre at z = 31500 mm.

Two field implementations are available:

### Uniform (default fallback)
If no `--field-map` argument is provided, the field volume contains a
uniform dipole field of **−0.15 T along X**.

### 3-D field map
If `--field-map <file>` is given, the simulation reads a tab-separated
text file and evaluates the field by trilinear interpolation.

The file format expected by `ShipFieldMap`:

```
X_SHiP   Y_SHiP   Z_SHiP   BX_SHiP   BY_SHiP   BZ_SHiP
-5000    -5000    -7500    0.001509  9.09e-05  0.000549
-5000    -5000    -7400    0.001508  9.72e-05  0.000543
...
```

- **Positions** (X, Y, Z) in mm, in the magnet's own coordinate system
  (origin at the magnet centre).
- **Field components** (BX, BY, BZ) in **Tesla**.
- A single header row with non-numeric column names is silently skipped.
- Lines starting with `#` are treated as comments.
- Must be a complete regular 3-D grid; the code infers (nX, nY, nZ) and
  the step sizes from the set of unique coordinates in the file.

The magnet-frame origin is mapped to **lab z = 31500 mm** (i.e. the centre
of the current field volume). Points outside the grid extent return B = 0.

---

## Project layout

```
StrawTracker/
├── CMakeLists.txt
├── straw_vis.mac
├── exec/
│   └── main.cpp                    Entry point and CLI parser
├── include/
│   ├── MaterialManager.h           Material definitions (incl. Aluminum)
│   ├── MagneticFieldRegion.h
│   ├── ShipFieldMap.h              G4MagneticField: 3-D map + interpolation
│   ├── StrawTrackerBuilder.h       GeoModel geometry builder
│   ├── TrackerDetectorConstruction.h
│   ├── TrackerHit.h                Hit data struct
│   ├── TrackerSD.h                 Sensitive detector
│   ├── TrackerRunAction.h
│   ├── TrackerEventAction.h
│   └── TrackerPrimaryGeneratorAction.h
└── src/
    ├── MaterialManager.cpp
    ├── MagneticFieldRegion.cpp
    ├── ShipFieldMap.cpp
    ├── StrawTrackerBuilder.cpp     Core geometry (GeoModel)
    ├── TrackerDetectorConstruction.cpp
    ├── TrackerSD.cpp
    ├── TrackerRunAction.cpp
    ├── TrackerEventAction.cpp
    └── TrackerPrimaryGeneratorAction.cpp
```

---

## Dependencies

| Package        | Minimum version | Notes                        |
|----------------|-----------------|------------------------------|
| CMake          | 3.16            |                              |
| GCC / Clang    | GCC 10 / Clang 11 | C++20 required             |
| Geant4         | 11.x            | Built with `ui_all vis_all`  |
| GeoModelCore   | 6.x             | Provides `GeoModelKernel`    |
| GeoModelG4     | 6.x             | Provides `GeoModel2G4`       |
| ROOT           | 6.x             | `Core`, `RIO`, `Tree`        |

---

## Building

```bash
mkdir build && cd build

cmake .. \
  -DCMAKE_PREFIX_PATH="/path/to/geant4;/path/to/geomodel;/path/to/root"

make -j$(nproc)
```

If the packages are already sourced via their setup scripts
(`thisgeant4.sh`, `thisroot.sh`, etc.), a plain `cmake ..` is sufficient.

---

## Running

### Visualisation only

```bash
./run_StrawTracker --visualize
```

### Simulate events (batch mode)

```bash
./run_StrawTracker --n-events 1000
```

### Export geometry to GDML

```bash
./run_StrawTracker --write-gdml
./run_StrawTracker --write-gdml --gdml-out my_geometry.gdml
```

### Run with a 3-D field map

```bash
./run_StrawTracker --n-events 1000 --field-map /path/to/field_map.txt
```

### Change the frame material

```bash
./run_StrawTracker --frame-material Aluminum   # default
./run_StrawTracker --frame-material Mylar
./run_StrawTracker --frame-material Air        # effectively no frame
```

### Full option reference

```
--n-events       <N>           Events to simulate            (default: 0)
--output         <file>        ROOT output file              (default: StrawTracker_hits.root)
--seed           <N>           Random seed (0 = auto)        (default: 0)
--particle       <n>           Geant4 particle name          (default: mu-)
--energy-MeV     <E>           Kinetic energy [MeV]          (default: 10000)
--pos-mm         <x> <y> <z>   Gun position [mm] lab frame   (default: 0 0 24000)
--dir            <x> <y> <z>   Direction unit vector         (default: 0 0 1)
--sigma-xy-mm    <s>           Gaussian beam spread X,Y [mm] (default: 0)
--field-map      <file>        Magnetic field-map text file
                               (empty = uniform -0.15 T fallback)
--frame-material <n>           Frame material (Aluminum,
                               Mylar, Air)                   (default: Aluminum)
--visualize                    Open interactive viewer
--vis-macro      <file>        Geant4 vis macro              (default: straw_vis.mac)
--write-gdml                   Export GDML
--gdml-out       <file>        GDML output file              (default: StrawTracker_geometry.gdml)
```

### Example commands

```bash
# 10 GeV muon, 1000 events, with a field map
./run_StrawTracker --n-events 1000 --particle mu- --energy-MeV 10000 \
    --field-map /data/magnet_field_map.txt

# Visualise with an aluminium frame
./run_StrawTracker --visualize --frame-material Aluminum

# Visualise and export GDML simultaneously
./run_StrawTracker --visualize --write-gdml

# Fixed seed for reproducible results
./run_StrawTracker --n-events 1000 --seed 42
```

---

## Output file structure

The ROOT file contains a TTree named `Events` with one row per event.
All branches are `std::vector<T>`; events with no hits have empty vectors.

| Branch       | Type            | Description                              |
|--------------|-----------------|------------------------------------------|
| `trackID`    | `vector<int>`   | Geant4 track ID                          |
| `stationID`  | `vector<int>`   | Station index (0–3)                      |
| `layerID`    | `vector<int>`   | Layer index within station (0–3)         |
| `subLayerID` | `vector<int>`   | Sub-layer index (0 = nominal, 1 = shifted) |
| `strawID`    | `vector<int>`   | Straw index within sub-layer (0–299)     |
| `edep`       | `vector<double>`| Energy deposit [MeV]                     |
| `x/y/z`      | `vector<double>`| Step centroid position [mm]              |
| `x/y/z_entry`| `vector<double>`| Track entry point into gas volume [mm]   |
| `x/y/z_exit` | `vector<double>`| Track exit point from gas volume [mm]    |

### Reading with uproot (Python)

```python
import uproot
import awkward as ak

with uproot.open("StrawTracker_hits.root") as f:
    tree = f["Events"]
    for batch in tree.iterate(step_size=100):
        print("station:", batch["stationID"])
        print("edep:",    batch["edep"], "MeV")
```

### Reading with ROOT macro (C++)

```cpp
void read_hits(const char* fname = "StrawTracker_hits.root") {
    TFile* f = TFile::Open(fname);
    TTree* t = f->Get<TTree>("Events");

    std::vector<int>*    stationID = nullptr;
    std::vector<double>* edep      = nullptr;
    t->SetBranchAddress("stationID", &stationID);
    t->SetBranchAddress("edep",      &edep);

    for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
        t->GetEntry(ev);
        std::cout << "Event " << ev << ": " << edep->size() << " hits\n";
    }
}
```

---

## Geometry notes

### Stereo angle implementation

The stereo rotation is applied to the **layer envelope volume** by the parent
station. Each layer logical volume is built axis-aligned (straws along X); the
`GeoTrf::RotateZ3D(±2.3° in rad)` transform is applied at placement time.
This means that within the layer's local frame, straws always run along X, which
simplifies drift-distance calculations. Because the frame is placed at z=0
inside the layer envelope, **the frame rotates with its view** — matching the
physical reality where the frame holds the straws.

### Staggering

The second sub-layer is displaced by `(+kStrawRadius, 0, +kStrawRadius)` in the
layer frame — i.e. half a pitch in Y and one radius in Z — giving the standard
double-staggered arrangement.

### Frame geometry

Each view's frame is built as a single GeoModel volume with shape
`GeoShapeSubtraction(outer_box, inner_hole_box)`, producing a hollow
rectangular frame. The frame material is set once per run via
`StrawTrackerBuilder::setFrameMaterial(name)` (called from
`TrackerDetectorConstruction` based on the `--frame-material` CLI argument).
Supported material names: `Aluminum` / `Aluminium` / `Al` / `Mylar` / `Air`.
Unknown names produce a warning and fall back to Aluminum.

### Wire (anode) geometry

The anode wire is not included in this geometry description. For simulations
requiring the wire, add a thin `GeoTube` (e.g. 25 µm tungsten) with `rMax ~
12.5e-3 mm` on the gas tube axis inside `buildStraw()`.

### Field map coordinate transform

The `ShipFieldMap` class treats the field-map origin as the magnet centre.
In the Geant4 world (which is centred at lab z = 31000 mm), the magnet
centre sits at z = +500 mm. The transform from a Geant4 world point to the
field-map frame is therefore `fm_point = g4_point - (0, 0, 500 mm)`. If
your magnet is positioned differently, change `MagneticFieldRegion::build`
where the `labZ` constant is set (currently `31500.0 * mm`).
