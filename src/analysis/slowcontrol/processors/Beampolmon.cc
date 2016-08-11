#include "Beampolmon.h"
#include <set>

using namespace std;
using namespace ant;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::processor;
using namespace ant::analysis::physics;

Processor::return_t Beampolmon::ProcessEventData(const TEventData& recon, manager_t& manager) {

    set<return_t> returnStates;

    returnStates.insert(Reference_1MHz.ProcessEventData(recon, manager));
    returnStates.insert(PbGlass.ProcessEventData(recon, manager));

    if (returnStates.size() != 1)
        throw runtime_error("Bug: Multiple return values for Scalar block.)");

    return *returnStates.begin();
}

void Beampolmon::PopQueue() {
    Reference_1MHz.PopQueue();
    PbGlass.PopQueue();
}

