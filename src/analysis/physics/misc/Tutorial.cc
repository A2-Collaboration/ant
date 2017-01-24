#include "Tutorial.h"

#include "base/Logger.h"

// use some namespaces (remember: only in implementation (aka .cc) files
// those statements are recommended to keep the following code not so namespace-clobbered
using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


Tutorial::Tutorial(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    LOG(INFO) << "We're such a cool Tutorial";

    // BinSettings nicely encapsulates this annoying (nbins,xlow,xhigh) stuff
    // and also provides generators/sanitization to counter binning artifacts
    // lets create 10 bins from 0 to 10
    BinSettings bins_nClusters(10,0,10);
    // shorter form: BinSettings bins_nClusters(10) (see constructors)

    // HistFac is a protected member of the base class "Physics"
    // use it to conveniently create histograms (and other ROOT objects) at the right location
    // using the make methods

    h_nClusters = HistFac.makeTH1D("Number of Clusters", // title
                                   "nClusters","#",      // xlabel, ylabel
                                   bins_nClusters,       // our binnings, may write directly BinSettings(10) here
                                   "h_nClusters"         // ROOT object name, auto-generated if omitted
                                   );
}

void Tutorial::ProcessEvent(const TEvent& event, manager_t&)
{
    h_nClusters->Fill(event.Reconstructed().Clusters.size());
}

// use the classes name to register the physics class inside Ant
// this is black macro magic what's used here...but it works :)
AUTO_REGISTER_PHYSICS(Tutorial)