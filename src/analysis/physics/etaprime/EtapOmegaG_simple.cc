#include "utils/ParticleTools.h"
#include "utils/Matcher.h"
#include "utils/Combinatorics.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "analysis/utils/uncertainties/Interpolated.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/std_ext/misc.h"

#include "physics/Physics.h"
#include "utils/ParticleTools.h"
#include "utils/fitter/TreeFitter.h"
#include "utils/MCWeighting.h"
#include "utils/A2GeoAcceptance.h"
#include "plot/PromptRandomHist.h"
#include "utils/TriggerSimulation.h"
#include "utils/ProtonPhotonCombs.h"

#include "tree/TSimpleParticle.h"
#include "base/ParticleTypeTree.h"
#include "base/WrapTTree.h"

class TH1D;
class TH2D;
class TH3D;

namespace ant {
namespace analysis {
namespace physics {

struct EtapOmegaG_simple : Physics {

    TH1D* h_Cuts = nullptr;

    TH1D* h_IM_Omega_true = nullptr;
    TH1D* h_IM_Etap_true = nullptr;

    TH1D* h_IM_2g;
    TH1D* h_IM_3g;
    TH1D* h_IM_4g;

    utils::TriggerSimulation triggersimu;
    PromptRandom::Switch promptrandom;

    const utils::UncertaintyModelPtr fitmodel_data;
    const utils::UncertaintyModelPtr fitmodel_mc;
    utils::KinFitter  kinfitter;

    EtapOmegaG_simple(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}}} // namespace ant::analysis::physics

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

EtapOmegaG_simple::EtapOmegaG_simple(const string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get()), // use setup for promptrandom windows
    fitmodel_data(// use Interpolated, based on Sergey's model
                  utils::UncertaintyModels::Interpolated::makeAndLoad(
                      utils::UncertaintyModels::Interpolated::Type_t::Data,
                      // use Sergey as starting point
                      make_shared<utils::UncertaintyModels::FitterSergey>()
                      )),
    fitmodel_mc(// use Interpolated, based on Sergey's model
                utils::UncertaintyModels::Interpolated::makeAndLoad(
                    utils::UncertaintyModels::Interpolated::Type_t::MC,
                    // use Sergey as starting point
                    make_shared<utils::UncertaintyModels::FitterSergey>()
                    )),
    kinfitter(nullptr, true)
{
    if(kinfitter.IsZVertexFitEnabled()) {
        const double z_vertex_sigma = 3.0;
        kinfitter.SetZVertexSigma(z_vertex_sigma);
        LOG(INFO) << "Fit Z vertex enabled with sigma=" << z_vertex_sigma;
    }

    h_Cuts = HistFac.makeTH1D("Cuts", "", "#", BinSettings(15),"h_Cuts");

    h_IM_Omega_true  = HistFac.makeTH1D("IM(3g) Omega True","IM / MeV","",{100,ParticleTypeDatabase::Omega.GetWindow(50)}, "h_IM_Omega_true");
    h_IM_Etap_true  = HistFac.makeTH1D("IM(4g) EtaPrime True","IM / MeV","",{100,ParticleTypeDatabase::EtaPrime.GetWindow(50)}, "h_IM_Etap_true");

    const AxisSettings axis_IM("IM / MeV",{500,0,1000});
    h_IM_2g = HistFac.makeTH1D("IM 2#gamma",axis_IM,"h_IM_2g");
    h_IM_3g = HistFac.makeTH1D("IM 3#gamma",axis_IM,"h_IM_3g");
    h_IM_4g = HistFac.makeTH1D("IM 4#gamma",axis_IM,"h_IM_4g");
}

void EtapOmegaG_simple::ProcessEvent(const TEvent& event, manager_t&)
{
    if(!triggersimu.ProcessEvent(event))
        h_Cuts->Fill("Triggersimu failed", 1.0);

    const TEventData& data = event.Reconstructed();

    const bool is_MC = data.ID.isSet(TID::Flags_t::MC);

    h_Cuts->Fill("Seen",1.0);

    // start now with some cuts
    if(!triggersimu.HasTriggered())
        return;
    h_Cuts->Fill("Triggered",1.0);

    if(data.Candidates.size()<3)
        return;
    h_Cuts->Fill("nCands>=3", 1.0);

    // etaprime physics always has something in TAPS
    // (the proton, by the way)
    // but remember: Backgrounds such as pi0pi0 have their proton in CB most likely
    // so it shouldn't be assumed that one of the TAPS clusters is the proton
    // (at least in the AntiPi0Pi0/AntiPi0Eta fits)
    bool haveTAPS = false;
    for(const auto& cand : data.Candidates.get_iter()) {
        if(cand->Detector & Detector_t::Type_t::TAPS) {
            haveTAPS = true;
        }
    }
    if(!haveTAPS)
        return;
    h_Cuts->Fill("1 in TAPS", 1.0);

    // this ensures the TParticlePtr (shared_ptr) are only made once
    // but do not allow photons with polar angle <7degree
    using particle_t = utils::ProtonPhotonCombs::comb_t;
    utils::ProtonPhotonCombs proton_photons(data.Candidates, [this] (particle_t& p) {
        auto it = p.Photons.begin();
        while(it != p.Photons.end()) {
            auto& photon = *it;
            if(std_ext::radian_to_degree(photon->Theta())<7) {
                it = p.Photons.erase(it);
            }
            else
                ++it;
        }
    });

    // set uncertainty model (maybe a bit ugly implemented here)
    kinfitter.SetUncertaintyModel(is_MC ? fitmodel_mc : fitmodel_data);


    for(const TTaggerHit& taggerhit : data.TaggerHits) {
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        auto combs = proton_photons()
                .Observe([this] (const std::string& s) { h_Cuts->Fill(s.c_str(), 1.0); })
                .FilterMult(4, 70.0)
                .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(350).Round());

        if(combs.empty())
            continue;

        auto kinFitProb = std_ext::NaN;
        vector<double> IM_2g(6, std_ext::NaN);
        vector<double> IM_3g(4, std_ext::NaN);
        vector<double> IM_4g(1, std_ext::NaN);

        for(auto& comb : combs) {

            auto result = kinfitter.DoFit(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

            if(result.Status != APLCON::Result_Status_t::Success)
                continue;

            if(!std_ext::copy_if_greater(kinFitProb, result.Probability))
                continue;

            const auto& photons = kinfitter.GetFittedPhotons();
            utils::ParticleTools::FillIMCombinations(IM_2g.begin(), 2, photons);
            utils::ParticleTools::FillIMCombinations(IM_3g.begin(), 3, photons);
            utils::ParticleTools::FillIMCombinations(IM_4g.begin(), 4, photons);
        }
    }
}

void EtapOmegaG_simple::ShowResult()
{
    canvas(GetName())
            << endc;
}

AUTO_REGISTER_PHYSICS(EtapOmegaG_simple)
