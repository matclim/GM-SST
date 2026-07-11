// ── Append this to src/StrawTrackerBuilder.cpp ───────────────────────────────
// Requires, near the top of the file:  #include <TFile.h>  #include <TTree.h>
//                                      #include <cmath>    #include <iostream>
//
// Writes a "Straws" TTree: one row per wire, in the WORLD frame (lab - worldZ),
// matching the frame the hits are recorded in. The reco reads this into
// ShipStrawGeometry to build one StrawSurface per wire.
//
// The placement mirrors buildStation/buildLayer/buildSubLayer exactly and
// reuses the same constants/helpers, so it stays in sync with the geometry.

void StrawTrackerBuilder::dumpStrawTable(const std::string& outFile,
                                         double worldZOriginMM) {
  const double layerGap   = 5.0;                              // mm (buildStation)
  const double layerPitch = 2.0 * layHalfZ_mm() + layerGap;   // 59 mm
  const double dzSub[2]   = { -(kStrawRadius + 0.55),         // SubLayer_0 nominal
                              +(kStrawRadius + 0.55) };       // SubLayer_1 shifted
  const double pitch      = 2.0 * kStrawRadius;               // 20 mm
  const double yStart     = -(kNStraws - 1) * 0.5 * pitch;    // -2990 mm
  const double halfLen    = 0.5 * kStrawLength;               // 2000 mm

  TFile f(outFile.c_str(), "RECREATE");
  TTree t("Straws", "Straw wire geometry (world frame, mm)");
  int    station = 0, layer = 0, subLayer = 0, straw = 0;
  double cx = 0, cy = 0, cz = 0, ux = 0, uy = 0, uz = 0;
  double radius = kStrawRadius, halfLength = halfLen;
  t.Branch("station",   &station);
  t.Branch("layer",     &layer);
  t.Branch("subLayer",  &subLayer);
  t.Branch("straw",     &straw);
  t.Branch("cx", &cx); t.Branch("cy", &cy); t.Branch("cz", &cz);
  t.Branch("ux", &ux); t.Branch("uy", &uy); t.Branch("uz", &uz);
  t.Branch("radius",     &radius);
  t.Branch("halfLength", &halfLength);

  for (int s = 0; s < 4; ++s) {
    const double zStation = kStationZ[s] - worldZOriginMM;
    for (int l = 0; l < kNLayers; ++l) {
      const double zLay  = -0.5 * (kNLayers - 1) * layerPitch + l * layerPitch;
      const double alpha = stereoSign(l) * kStereoAngle * M_PI / 180.0;  // rad
      const double ca = std::cos(alpha), sa = std::sin(alpha);
      for (int sub = 0; sub < 2; ++sub) {
        const double zWire    = zStation + zLay + dzSub[sub];
        const double yStagger = (sub == 1) ? kStrawRadius : 0.0;
        for (int i = 0; i < kNStraws; ++i) {
          const double yStraw = yStart + i * pitch + yStagger;
          station = s; layer = l; subLayer = sub; straw = i;
          // Wire centre = R_z(alpha) * (0, yStraw, 0) + (0,0,zWire).
          cx = -yStraw * sa;  cy = yStraw * ca;  cz = zWire;
          // Wire direction = R_z(alpha) * X  (straws run along X, then stereo).
          ux = ca;            uy = sa;           uz = 0.0;
          t.Fill();
        }
      }
    }
  }
  t.Write();
  std::cout << "[StrawTrackerBuilder] wrote " << t.GetEntries()
            << " straws to " << outFile << "\n";
  f.Close();
}
