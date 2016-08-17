#include "ExpTrigger.h"
#include <set>

using namespace std;
using namespace ant;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::processor;
using namespace ant::analysis::physics;

Processor::return_t ExpTrigger::ProcessEventData(const TEventData& recon, manager_t& manager) {

    set<return_t> returnStates;

    returnStates.insert(Reference_1MHz.ProcessEventData(recon, manager));
    returnStates.insert(LiveCounter.ProcessEventData(recon, manager));
    returnStates.insert(Trigger.ProcessEventData(recon,manager));
    returnStates.insert(L1Trigger.ProcessEventData(recon,manager));

    if (returnStates.size() != 1)
        throw runtime_error("Bug: Multiple return values for Scalar block.)");

    return *returnStates.begin();
}

void ExpTrigger::PopQueue() {
    Reference_1MHz.PopQueue();
    LiveCounter.PopQueue();
    Trigger.PopQueue();
    L1Trigger.PopQueue();
}

