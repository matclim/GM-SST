#pragma once
#include "G4VUserActionInitialization.hh"
#include <string>

class TrackerActionInitialization : public G4VUserActionInitialization {
public:
    TrackerActionInitialization(const std::string& outFile)
        : m_outFile(outFile) {}

    void BuildForMaster() const override;
    void Build() const override;

private:
    std::string m_outFile;
};
