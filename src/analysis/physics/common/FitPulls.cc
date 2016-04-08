#include "FitPulls.h"

#include "TH1D.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;



FitPulls::FitPulls(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{

}

void FitPulls::ProcessEvent(const TEvent& event, manager_t& manager)
{

}

void FitPulls::ShowResult()
{

}



AUTO_REGISTER_PHYSICS(FitPulls)