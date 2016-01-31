#include "DebugPhysics.h"

#include "base/Logger.h"
#include "expconfig/ExpConfig.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

DebugPhysics::DebugPhysics(const std::string& name, PhysOptPtr opts) :
    Physics(name, opts),
    writeEvents(opts->Get<bool>("WriteEvents", false)),
    requestSlowControl(opts->Get<bool>("RequestSlowControl", false))
{
}

DebugPhysics::~DebugPhysics() {}

void DebugPhysics::ProcessEvent(const TEvent& event, manager_t& manager)
{
    if(writeEvents) {
        manager.SaveEvent();
    }
    else {
        LOG(INFO) << event;
    }
}

void DebugPhysics::Finish()
{
    LOG(INFO) << "Nop";
}

void DebugPhysics::ShowResult()
{
    LOG(INFO) << "Nop";
}

void DebugPhysics::Initialize(input::SlowControl& slowcontrol)
{
    if(requestSlowControl)
        slowcontrol.FaradayCup.Request();
}



DebugPIDAlignment::DebugPIDAlignment(const std::string& name, PhysOptPtr opts):
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
