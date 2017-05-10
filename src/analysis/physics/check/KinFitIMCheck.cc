#include "KinFitIMCheck.h"

#include "expconfig/ExpConfig.h"
#include "utils/uncertainties/Interpolated.h"
#include "utils/ProtonPhotonCombs.h"
#include "utils/ParticleTools.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

KinFitIMCheck::KinFitIMCheck(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get()),
    fitmodel_data(// use Interpolated, based on Sergey's model
                  utils::UncertaintyModels::Interpolated::makeAndLoad(
                      utils::UncertaintyModels::Interpolated::Type_t::Data
                      )),
    fitmodel_mc(// use Interpolated, based on Sergey's model
                utils::UncertaintyModels::Interpolated::makeAndLoad(
                    utils::UncertaintyModels::Interpolated::Type_t::MC
                    )),
    fitter(nullptr, true)
{
    fitter.SetZVertexSigma(3);

    h_Steps = HistFac.makeTH1D("Steps","","",BinSettings(10),"h_Steps");
    t.CreateBranches(HistFac.makeTTree("t"));
}

void fill_IM(const TParticleList& photons,
             double& IM, vector<double>& IM_1, vector<double>& IM_2,vector<double>& IM_3,vector<double>& IM_4)
{
    LorentzVec photon_sum;
    for(auto& photon : photons)
        photon_sum += *photon;
    IM = photon_sum.M();

    const auto fill_IM_comb = [photons] (vector<double>& IM, int i) {

        const int nPhotons = photons.size();

        if(nPhotons-i>1) {
            IM.resize(std_ext::calcNchooseK(nPhotons, nPhotons-i));
            utils::ParticleTools::FillIMCombinations(IM.begin(), nPhotons-i, photons);
        }
        else
            IM.clear();
    };

    fill_IM_comb(IM_1, 1);
    fill_IM_comb(IM_2, 2);
    fill_IM_comb(IM_3, 3);
    fill_IM_comb(IM_4, 4);
}

void KinFitIMCheck::ProcessEvent(const TEvent& event, manager_t&)
{
    if(!triggersimu.ProcessEvent(event))
        h_Steps->Fill("Triggersimu failed", 1.0);

    if(!triggersimu.HasTriggered())
        return;

    const TEventData& data = event.Reconstructed();

    const bool is_MC = data.ID.isSet(TID::Flags_t::MC);
    fitter.SetUncertaintyModel(is_MC ? fitmodel_mc : fitmodel_data);


    if(data.Candidates.size() < 3)
        return;

    h_Steps->Fill("nCands>=3", 1.0);

    t.PIDSumE = 0;
    for(const TCluster& cl : data.Clusters) {
        if(cl.DetectorType == Detector_t::Type_t::PID) {
            t.PIDSumE += cl.Energy;
        }
    }

    t.nCB = 0;
    t.nTAPS = 0;
    t.CBVetoSumE = 0;

    for(const TCandidate& cand : data.Candidates) {
        if(cand.Detector & Detector_t::Type_t::CB) {
            t.nCB()++;
            t.CBVetoSumE() += cand.VetoEnergy;
        }
        else if(cand.Detector & Detector_t::Type_t::TAPS) {
            t.nTAPS()++;
        }
    }

    utils::ProtonPhotonCombs proton_photons(data.Candidates);

    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        h_Steps->Fill("Seen taggerhits",1.0);

        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;

        t.TaggW = promptrandom.FillWeight();
        t.TaggCh = taggerhit.Channel;

        auto filtered_combs = proton_photons()
                              .Observe([this] (const string& cut) { h_Steps->Fill(cut.c_str(), 1.0); }, "F ")
                              .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(350).Round());

        if(filtered_combs.empty()) {
            h_Steps->Fill("No combs left",1.0);
            continue;
        }

        // loop over the (filtered) proton combinations
        t.FitProb = std_ext::NaN;
        for(const auto& comb : filtered_combs) {

            const auto& result = fitter.DoFit(taggerhit.PhotonEnergy, comb.Proton, comb.Photons);

            if(result.Status != APLCON::Result_Status_t::Success)
                continue;
            if(!std_ext::copy_if_greater(t.FitProb, result.Probability))
                continue;

            t.ZVertex = fitter.GetFittedZVertex();

            fill_IM(comb.Photons,              t.IM_Raw, t.IM_Raw_1, t.IM_Raw_2, t.IM_Raw_3, t.IM_Raw_4);
            fill_IM(fitter.GetFittedPhotons(), t.IM,     t.IM_1,     t.IM_2,     t.IM_3,     t.IM_4);
        }


        if(t.FitProb > 0.001) {
            h_Steps->Fill("TreeFill",1.0);
            t.Tree->Fill();
        }
    }
}

void KinFitIMCheck::ShowResult()
{
    canvas(GetName())
            << h_Steps
            << endc;
}


AUTO_REGISTER_PHYSICS(KinFitIMCheck)