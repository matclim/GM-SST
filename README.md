# ship_reco — standalone ACTS tracking & vertexing

Reconstruction library that reads the Geant4 `Events` tree, fits muon tracks
through the spectrometer field map with the ACTS KalmanFitter, and (later)
reconstructs decay vertices.

## Layout
    field/       ROOT field map -> Acts::MagneticFieldProvider          [ready]
    geometry/    straw layers -> free PlaneSurfaces (DirectNavigator)   [fill numbers]
    edm/         SourceLink, 1-D measurement, calibrator, accessor      [API-sensitive]
    io/          Events tree -> in-memory hits                          [ready]
    seeding/     upstream-station straight-line seed + p guess          [ready]
    fitting/     KalmanFitter over the surface sequence                 [API-sensitive]
    vertexing/   Billoir / AMVF on fitted tracks                        [stub, do last]
    app/         run_reco driver

## Build
Configure your project pointing at your ACTS install, e.g.
    cmake -B build -S . -DActs_DIR=$HOME/Software/acts/install/lib/cmake/Acts
    cmake --build build -j

## Run
    ./build/run_reco --hits field_muon.root \
                     --field ../../fieldmaps/NEWMAPS/2026_07_02_MainSpectrometerField_V21_2455.root \
                     --p-guess 10000 --field-origin-z 500

## Bring-up order (each step fails independently)
1. Compile `field` only, confirm it links (already validated numerically).
2. Add `geometry` + `io` + `seeding`; run with a truth-position measurement and
   check the propagator visits the right surfaces.
3. Add `edm` + `fitting`; get a KalmanFitter result; compare <p> to the gun energy.
4. Replace plane/position measurements with real straw drift radii + L/R ambiguity,
   and switch the finder to the CombinatorialKalmanFilter.
5. Enable `vertexing` for multi-track events.

## Notes on the two frames
`--field-origin-z 500` and `defaultLayers(31000)` both assume the ACTS global
frame == the Geant4 world (centre lab z = 31000 mm). Keep them consistent.

## Lines marked `// ACTS-API`
These are the calls most likely to differ between ACTS releases (measurement
calibration, KalmanFitter extensions/options, vertexing). If one doesn't
resolve, mirror the equivalent from your ACTS `Examples/` sources.
