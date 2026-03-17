#pragma once
// MaterialManager.h
// Defines and caches all materials used by the straw tracker.

#include "GeoModelKernel/GeoMaterial.h"
#include "GeoModelKernel/GeoElement.h"

class MaterialManager {
public:
    static MaterialManager& instance();

    // Access individual materials (created once, returned by pointer).
    GeoMaterial* ArCO2();     // 70% Ar + 30% CO2 by mass
    GeoMaterial* Mylar();     // Polyethylene terephthalate (C10H8O4)
    GeoMaterial* Air();       // Dry air (world fill)

private:
    MaterialManager();

    GeoMaterial* m_ArCO2  {nullptr};
    GeoMaterial* m_Mylar  {nullptr};
    GeoMaterial* m_Air    {nullptr};
};
