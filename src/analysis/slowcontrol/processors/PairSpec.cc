#include "PairSpec.h"

#include <set>

using namespace std;
using namespace ant;
using namespace ant::analysis::slowcontrol;
using namespace ant::analysis::slowcontrol::processor;
using namespace ant::analysis::physics;


Processor::return_t Beam::ProcessEventData(const TEventData& recon, manager_t& manager)
{
     set<return_t> returnStates;

     returnStates.insert(IonChamber.ProcessEventData(recon, manager));
     returnStates.insert(PairSpecGate.ProcessEventData(recon, manager));

     if (returnStates.size() != 1)
         throw runtime_error("Bug: Multiple return values for Scalar block.)");

     return *returnStates.begin();
}

void Beam::PopQueue()
{
    IonChamber.PopQueue();
    PairSpecGate.PopQueue();
}
