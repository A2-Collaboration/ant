#include "DebugPhysics.h"

#include "slowcontrol/SlowControlVariables.h"

#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

DebugPhysics::DebugPhysics(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    writeEvents(opts->Get<unsigned>("WriteEvents", 0)),
    keepReadHits(opts->Get<bool>("KeepReadHits", false)),
    requestSlowControl(opts->Get<bool>("RequestSlowControl", false))
{
    slowcontrol::Variables::TaggerScalers->Request();
}

DebugPhysics::~DebugPhysics() {}

void DebugPhysics::ProcessEvent(const TEvent& event, manager_t& manager)
{
    if(writeEvents>0 && seenEvents % writeEvents == 0) {
        manager.SaveEvent();
        if(keepReadHits)
            manager.KeepDetectorReadHits();
    }
    else if(!writeEvents) {
        LOG(INFO) << event;
    }
    seenEvents++;
}

void DebugPhysics::Finish()
{
    LOG(INFO) << "Seen " << seenEvents << " Events";
}

void DebugPhysics::ShowResult()
{
    LOG(INFO) << "Nop";
}

void DebugPhysics::Initialize(slowcontrol::SlowControl& slowcontrol)
{
    if(requestSlowControl)
        slowcontrol.FaradayCup.Request();
}



DebugPIDAlignment::DebugPIDAlignment(const std::string& name, OptionsPtr opts):
    Physics(name, opts)
{
    const BinSettings bins(360,-180,180);
    angles = HistFac.makeTH2D("PID Phi angles","MCTrue #phi","Rec #phi",bins,bins);
}

DebugPIDAlignment::~DebugPIDAlignment()
{

}

void DebugPIDAlignment::ProcessEvent(const TEvent& event, manager_t&)
{
    if(event.MCTrue->Particles.GetAll().size() == 1) {
        const auto mctrue_phi = event.MCTrue->Particles.GetAll().front()->Phi() * TMath::RadToDeg();

        for(const TCandidatePtr& cand : event.Reconstructed->Candidates) {
            for(const TClusterPtr& c : cand->Clusters) {
                if(c->DetectorType == Detector_t::Type_t::PID) {
                    angles->Fill(mctrue_phi, c->Position.Phi()* TMath::RadToDeg());
                }
            }
        }

        for(const TClusterPtr& c : event.Reconstructed->Clusters) {
            if(c->DetectorType == Detector_t::Type_t::PID) {
                angles->Fill(mctrue_phi, c->Position.Phi()* TMath::RadToDeg());
            }
        }
    }
}

void DebugPIDAlignment::ShowResult()
{
    canvas("PID angles") << drawoption("colz") << angles << endc;

}

AUTO_REGISTER_PHYSICS(DebugPhysics)
AUTO_REGISTER_PHYSICS(DebugPIDAlignment)
