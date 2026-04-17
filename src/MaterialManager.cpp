// MaterialManager.cpp
#include "MaterialManager.h"

#include "GeoModelKernel/GeoMaterial.h"
#include "GeoModelKernel/GeoElement.h"

#include <algorithm>
#include <cctype>
#include <iostream>

namespace GU = GeoModelKernelUnits;

MaterialManager& MaterialManager::instance() {
    static MaterialManager inst;
    return inst;
}

MaterialManager::MaterialManager() {}

// ── Ar/CO2 70/30 by mass ──────────────────────────────────────────────────────
GeoMaterial* MaterialManager::ArCO2() {
    if (m_ArCO2) return m_ArCO2;

    // Elements
    auto* Ar = new GeoElement("Argon",   "Ar", 18, 39.948 * GU::g / GU::mole);
    auto* C  = new GeoElement("Carbon",  "C",   6, 12.011 * GU::g / GU::mole);
    auto* O  = new GeoElement("Oxygen",  "O",   8, 15.999 * GU::g / GU::mole);

    // Density of Ar/CO2 70/30 at STP ≈ 1.56 kg/m³ → 1.56e-3 g/cm³
    const double density = 1.56e-3 * GU::g / GU::cm3;

    m_ArCO2 = new GeoMaterial("ArCO2_70_30", density);

    const double mCO2     = 44.009;
    const double frac_CO2 = 0.30;
    const double frac_Ar  = 0.70;

    m_ArCO2->add(Ar, frac_Ar);
    m_ArCO2->add(C,  frac_CO2 * 12.011 / mCO2);
    m_ArCO2->add(O,  frac_CO2 * 2.0 * 15.999 / mCO2);

    m_ArCO2->lock();
    return m_ArCO2;
}

// ── Mylar (PET, C10H8O4) ──────────────────────────────────────────────────────
GeoMaterial* MaterialManager::Mylar() {
    if (m_Mylar) return m_Mylar;

    auto* H = new GeoElement("Hydrogen", "H",  1,  1.008  * GU::g / GU::mole);
    auto* C = new GeoElement("Carbon",   "C",  6, 12.011  * GU::g / GU::mole);
    auto* O = new GeoElement("Oxygen",   "O",  8, 15.999  * GU::g / GU::mole);

    const double density = 1.39 * GU::g / GU::cm3;
    const double mPET    = 10*12.011 + 8*1.008 + 4*15.999;

    m_Mylar = new GeoMaterial("Mylar", density);
    m_Mylar->add(C, 10*12.011 / mPET);
    m_Mylar->add(H,  8* 1.008 / mPET);
    m_Mylar->add(O,  4*15.999 / mPET);
    m_Mylar->lock();
    return m_Mylar;
}

// ── Air ───────────────────────────────────────────────────────────────────────
GeoMaterial* MaterialManager::Air() {
    if (m_Air) return m_Air;

    auto* N = new GeoElement("Nitrogen", "N",  7, 14.007 * GU::g / GU::mole);
    auto* O = new GeoElement("Oxygen",   "O",  8, 15.999 * GU::g / GU::mole);

    const double density = 1.205e-3 * GU::g / GU::cm3;

    m_Air = new GeoMaterial("Air", density);
    m_Air->add(N, 0.755);
    m_Air->add(O, 0.232);
    m_Air->lock();
    return m_Air;
}

// ── Aluminum (structural frames) ──────────────────────────────────────────────
GeoMaterial* MaterialManager::Aluminum() {
    if (m_Aluminum) return m_Aluminum;

    auto* Al = new GeoElement("Aluminium", "Al", 13, 26.9815 * GU::g / GU::mole);

    const double density = 2.70 * GU::g / GU::cm3;   // pure Al

    m_Aluminum = new GeoMaterial("Aluminum", density);
    m_Aluminum->add(Al, 1.0);
    m_Aluminum->lock();
    return m_Aluminum;
}

// ── Lookup by name (case-insensitive) ────────────────────────────────────────
GeoMaterial* MaterialManager::frameMaterialByName(const std::string& name) {
    std::string key = name;
    std::transform(key.begin(), key.end(), key.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (key == "al" || key == "aluminum" || key == "aluminium") return Aluminum();
    if (key == "mylar" || key == "pet")                         return Mylar();
    if (key == "air")                                           return Air();

    std::cerr << "[MaterialManager] Unknown frame material '" << name
              << "', falling back to Aluminum.\n";
    return Aluminum();
}
