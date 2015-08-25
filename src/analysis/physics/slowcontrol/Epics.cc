#include "Epics.h"

using namespace std;
using namespace ant::analysis::slowcontrol;

EpicsPV::MyGetter::MyGetter(const string& pvname):pv(pvname) {}

void EpicsPV::MyGetter::Process(double d) {v=d;}

std::list<ant::TSlowControl::Key> EpicsPV::MyGetter::GetRequiredKeys() const {
    return { TSlowControl::Key(TSlowControl::Type_t::EpicsScaler, pv) };
}


TotalLivetime::TotalLivetime(Receiver* receiver): EpicsPV(receiver,"TRIG:TotalLivetime") {}


EpicsPV::EpicsPV(Receiver* receiver, const string& pvname):Variable(receiver), pv(pvname) {}
