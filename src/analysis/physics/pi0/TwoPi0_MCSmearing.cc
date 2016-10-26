#include "TwoPi0_MCSmearing.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"
#include "utils/Uncertainties.h"
#include "base/Logger.h"
#include "expconfig/ExpConfig.h"
#include "tree/TParticle.h"
#include "tree/TCandidate.h"
#include "tree/TCluster.h"

#include "TH1D.h"
#include "TTree.h"

#include <memory>
#include <cassert>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;



TwoPi0_MCSmearing::TwoPi0_MCSmearing(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{

}

TwoPi0_MCSmearing::~TwoPi0_MCSmearing()
{}

void TwoPi0_MCSmearing::ProcessEvent(const TEvent& event, manager_t&)
{

}


AUTO_REGISTER_PHYSICS(TwoPi0_MCSmearing)
