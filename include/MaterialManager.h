#pragma once
// MaterialManager.h
// Defines and caches all materials used by the straw tracker.

#include "GeoModelKernel/GeoMaterial.h"
#include "GeoModelKernel/GeoElement.h"

#include <string>

class MaterialManager {
public:
    static MaterialManager& instance();

    // Access individual materials (created once, returned by pointer).
    GeoMaterial* ArCO2();     // 70% Ar + 30% CO2 by mass
    GeoMaterial* Mylar();     // Polyethylene terephthalate (C10H8O4)
    GeoMaterial* Air();       // Dry air (world fill)
    GeoMaterial* Aluminum();  // Structural aluminium (frames)

    // Lookup a frame material by name. Case-insensitive.
    // Falls back to Aluminum() with a warning if the name is unknown.
    // Supported names: "Al", "Aluminum", "Aluminium", "Mylar", "Air".
    GeoMaterial* frameMaterialByName(const std::string& name);

private:
    MaterialManager();

    GeoMaterial* m_ArCO2     {nullptr};
    GeoMaterial* m_Mylar     {nullptr};
    GeoMaterial* m_Air       {nullptr};
    GeoMaterial* m_Aluminum  {nullptr};
};
