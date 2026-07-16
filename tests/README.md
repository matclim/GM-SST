# tests/

Analysis and end-to-end scripts for the StrawTracker + ACTS reco chain.
Assumes the project is built in `../build` (run_StrawTracker, run_reco).

## End-to-end
- **run_all.sh** — generates 3 muon samples in slightly different directions
  (straight / tiltx / tilty), 200 events each, through the field map, then
  reconstructs each into `reco_<name>.root`. Resolves paths relative to itself,
  so it runs from anywhere:

      ./tests/run_all.sh
      FIELD=/path/to/other_map.root ./tests/run_all.sh   # override the map
      OUTDIR=/tmp/recotest ./tests/run_all.sh             # override output dir

## Field-map checks (ROOT macros)
- **check_orientation.C** — verifies a field map's axes, units and dominant
  component:  `root -l 'check_orientation.C("map.root", 0.15)'`
- **check_bending.C** — measures the muon bend angle from a reco/hits file and
  compares to the analytic field integral:
  `root -l 'check_bending.C("field_muon.root","map.root",10000)'`
- **plot_field_vs_z.C** — plots Bx,By vs z on the axis:
  `root -l 'plot_field_vs_z.C("map.root")'`

## Looking at reco output
    root -l reco_straight.root
    root [1] Tracks->Draw("pRes")             // (p_fit - p_true)/p_true
    root [2] Tracks->Draw("p")                // fitted |p| [GeV]
    root [3] Tracks->Draw("chi2/(nMeas-5)")   // chi2 / dof
