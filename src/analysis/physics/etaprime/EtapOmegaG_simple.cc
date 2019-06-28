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

    TH1D* h_Cuts;

    TH1D* h_IM_Omega_true;
    TH1D* h_IM_Etap_true;

    TH1D* h_KinFitProb;
    TH2D* h_ZVertex;

    TH1D* h_IM_2g;
    TH1D* h_IM_3g;
    TH1D* h_IM_4g;

    TH1D* h_IM_Pi0;
    TH1D* h_IM_Pi0g_max;

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
//    fitmodel_data(make_shared<utils::UncertaintyModels::FitterSergey>()),
//    fitmodel_mc(make_shared<utils::UncertaintyModels::FitterSergey>()),
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

    h_KinFitProb = HistFac.makeTH1D("KinFitProb",{"p",{100,{0,1}}},"h_KinFitProb");

    const BinSettings bins_zvertex(50,-10,10);
    h_ZVertex = HistFac.makeTH2D("ZVertex Fitted vs. True",
                                 {"z_{true} / cm", bins_zvertex},
                                 {"z_{fitted} / cm", bins_zvertex},
                                 "h_ZVertex");

    h_IM_2g = HistFac.makeTH1D("IM All 2#gamma combinations",{"IM(2#gamma) / MeV",{200,ParticleTypeDatabase::Pi0.GetWindow(200)}},"h_IM_2g");
    h_IM_3g = HistFac.makeTH1D("IM All 3#gamma combinations",{"IM(3#gamma) / MeV",{200,ParticleTypeDatabase::Omega.GetWindow(200)}},"h_IM_3g");
    h_IM_4g = HistFac.makeTH1D("IM 4#gamma",{"IM(4#gamma) / MeV",{200,ParticleTypeDatabase::EtaPrime.GetWindow(200)}},"h_IM_4g");

    h_IM_Pi0 = HistFac.makeTH1D("IM Best #pi^{0}",{"IM(#pi^{0}}) / MeV",{200,ParticleTypeDatabase::Pi0.GetWindow(200)}},"h_IM_Pi0");
    h_IM_Pi0g_max = HistFac.makeTH1D("IM Max #pi^{0}g",{"IM(#pi^{0}g)_{max} / MeV",{200,ParticleTypeDatabase::Omega.GetWindow(200)}},"h_IM_Pi0g_max");
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

    if(data.Candidates.size()<5)
        return;
    h_Cuts->Fill("nCands>=5", 1.0);

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
                .Observe([this] (const std::string& s) { h_Cuts->Fill(s.c_str(), 1.0); },"F ")
                .FilterMult(4, 0.0)
                .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(350).Round());

        if(combs.empty())
            continue;

        auto kinFitProb = std_ext::NaN;
        auto best_zvertex = std_ext::NaN;
        TParticleList best_photons;

        for(auto& comb : combs) {
            auto result = kinfitter.DoFit(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);
            if(result.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(kinFitProb, result.Probability))
                continue;
            best_photons = kinfitter.GetFittedPhotons();
            best_zvertex = kinfitter.GetFittedZVertex();
        }

        if(!(kinFitProb>0.01))
            continue;


        h_KinFitProb->Fill(kinFitProb, promptrandom.FillWeight());
        h_ZVertex->Fill(event.MCTrue().Target.Vertex.z, best_zvertex, promptrandom.FillWeight());

        h_Cuts->Fill("KinFitProb>0.01",1.0);

        // fill stupid combinatorics
        utils::ParticleTools::FillIMCombinations(h_IM_2g, 2, best_photons, promptrandom.FillWeight());
        utils::ParticleTools::FillIMCombinations(h_IM_3g, 3, best_photons, promptrandom.FillWeight());
        utils::ParticleTools::FillIMCombinations(h_IM_4g, 4, best_photons, promptrandom.FillWeight());

        // find best photons for pi0 by looking at IM difference (basically chi2 test)
        TParticleList pi0_photons(2);
        {
            auto pi0_IM_diff = std_ext::NaN;
            for( auto comb = utils::makeCombination(best_photons,2); !comb.done(); ++comb) {
                const auto pi0 = *comb.at(0) + *comb.at(1);
                const auto IM_diff = std::fabs(pi0.M() - ParticleTypeDatabase::Pi0.Mass());
                if(!std_ext::copy_if_better(pi0_IM_diff, IM_diff, std::less<double>()))
                    continue;
                pi0_photons.front() = comb.at(0);
                pi0_photons.back() = comb.at(1);
            }
        }

        const auto pi0 = *pi0_photons.front() + *pi0_photons.back();
        h_IM_Pi0->Fill(pi0.M(), promptrandom.FillWeight());

        // find the two remaining bachelor photons (not beloning to the pi0)
        TParticleList bachelor_photons;
        for(const auto& p : best_photons) {
            if(!std_ext::contains(pi0_photons, p)) {
                bachelor_photons.push_back(p);
            }
        }

        // build the two possible pi0g IMs
        const auto IM_Pi0g_1 = (pi0 + *bachelor_photons.front()).M();
        const auto IM_Pi0g_2 = (pi0 + *bachelor_photons.back()).M();
        // higher one expected to be omega
        h_IM_Pi0g_max->Fill(IM_Pi0g_1>IM_Pi0g_2?IM_Pi0g_1:IM_Pi0g_2, promptrandom.FillWeight());
    }
}

void EtapOmegaG_simple::ShowResult()
{
    canvas(GetName())
            << h_Cuts << h_KinFitProb << drawoption("colz") << h_ZVertex
            << endr
            << h_IM_2g << h_IM_3g << h_IM_4g
            << endr
            << h_IM_Pi0 << h_IM_Pi0g_max
            << endc;
}

AUTO_REGISTER_PHYSICS(EtapOmegaG_simple)
