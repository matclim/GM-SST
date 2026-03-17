#include "TrackerActionInitialization.h"
#include "TrackerRunAction.h"
#include "TrackerEventAction.h"
#include "TrackerPrimaryGeneratorAction.h"

void TrackerActionInitialization::BuildForMaster() const {
    // Master thread only needs the run action
    SetUserAction(new TrackerRunAction(nullptr, m_outFile));
}

void TrackerActionInitialization::Build() const {
    auto* evtAction = new TrackerEventAction();
    SetUserAction(new TrackerPrimaryGeneratorAction());
    SetUserAction(evtAction);
    SetUserAction(new TrackerRunAction(evtAction, m_outFile));
}
