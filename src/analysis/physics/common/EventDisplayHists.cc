#include "EventDisplayHists.h"

#include "utils/ParticleTools.h"

#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"

#include "root-addons/cbtaps_display/TH2TAPS.h"
#include "root-addons/cbtaps_display/TH2CB.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/TAPS.h"

#include "TCanvas.h"
#include "TMarker.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis::physics;
using namespace std;


EventDisplayHists::EventDisplayHists(const string& name, OptionsPtr opts):
    Physics(name, opts)
{
    const auto taps = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
    tapsZ = taps->GetZPosition();
}

EventDisplayHists::~EventDisplayHists()
{}

void EventDisplayHists::ProcessEvent(const TEvent& event, manager_t&)
{

    const auto& candidates = event.Reconstructed().Candidates;

    auto taps_cands = candidates.get_ptr_list(
                          [] (const TCandidate& cand) {
        return cand.Detector & Detector_t::Type_t::TAPS;
    });

    if(taps_cands.empty())
        return;

    auto tapsCal = HistFac.make<TH2TAPS>(formatter() << "taps" << n++, formatter() << taps_cands.size() << " Candidates (TAPS)");

    for(const auto& c : taps_cands) {
        const auto cluster = c->FindCaloCluster();

        for(const auto& hit : cluster->Hits) {
            tapsCal->SetElement(hit.Channel, hit.Energy);
        }

        tapsCal->CreateMarker(cluster->Position.XY(), kFullCircle, kOpenCircle);

        tapsCal->CreateMarker(cluster->CentralElement);

    }
    auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);

    const auto true_particles = mctrue_particles.GetAll();

    for(const auto& p : true_particles) {
        if(p->Theta() < degree_to_radian(30.0)) {
            const auto xy = p->p.XY() * tapsZ / p->p.z;

            tapsCal->CreateMarker(xy, kFullSquare, kOpenSquare);

        }
    }
}


AUTO_REGISTER_PHYSICS(EventDisplayHists)
