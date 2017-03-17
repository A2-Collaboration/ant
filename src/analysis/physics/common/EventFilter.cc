#include "EventFilter.h"

#include "TH1D.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "utils/ProtonPermutation.h"


using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;



void SetBinLabel(TH1* hist, const int bin, const std::string& label) {
    hist->GetXaxis()->SetBinLabel(bin, label.c_str());
}

LorentzVec sum_particles(const TParticleList& particles) {

    LorentzVec lv;

    for(const auto& p : particles) {
        lv += *p;
    }

    return lv;
}

bool EventFilter::checkCoplanarity(const TCandidateList &cands, const double maxangle)
{

    for(utils::ProtonPermutation perm(cands.get_ptr_list()); perm.Good(); perm.Next()) {
        const LorentzVec photons = sum_particles(perm.Photons());
        const auto& proton = perm.Proton();

        if(fabs(vec2::Phi_mpi_pi(M_PI + proton->Phi() - photons.Phi())) < maxangle)
            return true;
    }

    return false;

}

EventFilter::EventFilter(const string& name, OptionsPtr opts):
    Physics(name, opts),

    nCands(opts->Get<interval<size_t>>("nCands", {0, 100})),
    CBEsum(opts->Get<double>("CBEsum",    550.0)),
    maxCoplAngle(degree_to_radian(opts->Get<double>("CoplAngle", 15.0)))

{
    steps = HistFac.makeTH1D("Filter Steps", "", "", BinSettings(4));
    SetBinLabel(steps,1, formatter() << "Input");
    SetBinLabel(steps,2, formatter() << "CBEsum > " << CBEsum);
    SetBinLabel(steps,3, formatter() << "nCands in " << nCands);
    SetBinLabel(steps,4, formatter() << "Copl angle < " << radian_to_degree(maxCoplAngle) << "#circ");

    VLOG(3) << "CBEsum > " << CBEsum;
    VLOG(3) << "nCands in " << nCands;
    VLOG(3) << "Coplanarity angle < " << radian_to_degree(maxCoplAngle) << " degree";
}

EventFilter::~EventFilter()
{}



void EventFilter::ProcessEvent(const TEvent& event, manager_t& manager)
{
    triggersimu.ProcessEvent(event);

    const auto& data = event.Reconstructed();

    steps->SetBinContent(1, steps->GetBinContent(1)+1);

    if(triggersimu.GetCBEnergySum() < CBEsum)
        return;

    steps->SetBinContent(2, steps->GetBinContent(2)+1);

    if(!nCands.Contains(data.Candidates.size()))
        return;

    steps->SetBinContent(3, steps->GetBinContent(3)+1);

    if(!checkCoplanarity(data.Candidates, maxCoplAngle))
        return;

    steps->SetBinContent(4, steps->GetBinContent(4)+1);

    manager.SaveEvent();
}

AUTO_REGISTER_PHYSICS(EventFilter)
