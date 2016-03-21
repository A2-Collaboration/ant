#include "KinFitPi0.h"

#include "utils/particle_tools.h"
#include "utils/matcher.h"

#include "expconfig/ExpConfig.h"

#include "TH1D.h"
#include "TTree.h"

#include <memory>
#include <cassert>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;



KinFitPi0::KinFitPi0(const string& name, OptionsPtr opts) :
    Physics(name, opts)
{
    for(unsigned mult=1;mult<=opts->Get<unsigned>("nPi0",1);mult++) {
        multiPi0.emplace_back(std_ext::make_unique<MultiPi0>(HistFac, mult));
    }
}

void KinFitPi0::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& data = event.Reconstructed();
    for(auto& m : multiPi0)
        m->ProcessData(data);
}

void KinFitPi0::ShowResult()
{
    for(auto& m : multiPi0)
        m->ShowResult();
}




KinFitPi0::MultiPi0::MultiPi0(HistogramFactory& histFac, unsigned nPi0) :
    multiplicity(nPi0),
    fitter(std_ext::formatter() << multiplicity << "Pi0", 2*multiplicity),
    h_missingmass(promptrandom),
    h_fitprobability(promptrandom),
    IM_2g_byFit(promptrandom),
    IM_2g_fitted(promptrandom)
{
    std::string multiplicity_str = std_ext::formatter() << "m" << multiplicity << "Pi0";
    HistogramFactory HistFac(multiplicity_str, histFac, multiplicity_str);

    promptrandom.AddPromptRange({-2.5,1.5}); // slight offset due to CBAvgTime reference
    promptrandom.AddRandomRange({-50,-10});  // just ensure to be way off prompt peak
    promptrandom.AddRandomRange({  10,50});

    steps = HistFac.makeTH1D("Steps","","#",BinSettings(15),"steps");

    Proton_Coplanarity = HistFac.makeTH1D("p Coplanarity","#delta#phi / degree","",BinSettings(400,-180,180),"Proton_Coplanarity");

    h_taggtime = HistFac.makeTH1D("Tagged Time","t / ns", "", BinSettings(300,-60,60), "h_taggtime");
    h_missingmass.MakeHistograms(HistFac, "h_missingmass","Missing Mass",BinSettings(400,400, 1400),"MM / MeV","#");
    h_fitprobability.MakeHistograms(HistFac, "fit_probability","KinFitter probability",BinSettings(150,0,1),"p","#");

    BinSettings bins_IM(1400,0,1400);

    IM_2g_byFit.MakeHistograms(HistFac, "IM_2g_byFit","IM 2#gamma by Fit",bins_IM,"IM / MeV","#");
    IM_2g_fitted.MakeHistograms(HistFac, "IM_2g_fitted","IM 2#gamma fitted",bins_IM,"IM / MeV","#");


    const auto setup = ant::ExpConfig::Setup::GetLastFound();
    fitter.LoadSigmaData(setup->GetPhysicsFilesDirectory()+"/FitterSigmas.root");
}

void KinFitPi0::MultiPi0::ProcessData(const TEventData& data)
{
    const auto nPhotons_expected = multiplicity*2;

    steps->Fill("Seen",1);

    const auto& cands = data.Candidates;
    const auto nCandidates = cands.size();
    const auto nCandidates_expected = nPhotons_expected+1;
    if(nCandidates != nCandidates_expected)
        return;
    std::string nCandidates_cutstr = std_ext::formatter() << "nCandidates==" << nCandidates_expected;
    steps->Fill(nCandidates_cutstr.c_str(),1);

    // use any candidate as proton, and do the analysis (ignore ParticleID stuff)

    for(auto i_proton : cands.get_iter()) {

        const auto proton = std::make_shared<TParticle>(ParticleTypeDatabase::Proton, i_proton);
        std::vector<TParticlePtr> photons;
        unsigned nPhotons_CB = 0;
        unsigned nPhotons_TAPS = 0;
        bool touchesHole = false;
        for(auto i_photon : cands.get_iter()) {
            if(i_photon == i_proton)
                continue;
            photons.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon, i_photon));

            if(i_photon->Detector & Detector_t::Type_t::CB)
                nPhotons_CB++;
            else if(i_photon->Detector & Detector_t::Type_t::TAPS)
                nPhotons_TAPS++;

            if(i_photon->FindCaloCluster()->HasFlag(TCluster::Flags_t::TouchesHole))
                touchesHole = true;
        }
        assert(photons.size() == nPhotons_expected);

        if(nPhotons_CB != nPhotons_expected) {
            steps->Fill("Photons in TAPS", 1.0);
            continue;
        }

        if(touchesHole) {
            steps->Fill("TouchesHole", 1.0);
            continue;
        }

        LorentzVec photon_sum(0,0,0,0);
        for(const auto& p : photons) {
            photon_sum += *p;
        }

        // proton coplanarity check

        const double d_phi = std_ext::radian_to_degree(vec2::Phi_mpi_pi(proton->Phi()-photon_sum.Phi() - M_PI ));
        Proton_Coplanarity->Fill(d_phi);

        // iterate over tagger hits

        for(const TTaggerHit& taggerhit : data.TaggerHits) {

            steps->Fill("Seen taggerhits",1.0);
            const auto taggtime = taggerhit.Time - data.Trigger.CBTiming;
            promptrandom.SetTaggerHit(taggtime);
            if(promptrandom.State() == PromptRandom::Case::Outside)
                continue;
            h_taggtime->Fill(taggtime);

            // simple missing mass cut
            const LorentzVec beam_target = taggerhit.GetPhotonBeam() + LorentzVec(0, 0, 0, ParticleTypeDatabase::Proton.Mass());
            const LorentzVec missing = beam_target - photon_sum;
            const double missing_mass = missing.M();

            h_missingmass.Fill(missing_mass);
            const interval<double> MM_cut(800, 1050);
            if(!MM_cut.Contains(missing_mass))
                continue;

            steps->Fill("Fitted taggerhits",1.0);

            // more sophisticated fitter
            fitter.SetEgammaBeam(taggerhit.PhotonEnergy);
            fitter.SetProton(proton);
            fitter.SetPhotons(photons);
            auto fit_result = fitter.DoFit();

            if(fit_result.Status == APLCON::Result_Status_t::Success) {
                steps->Fill("Fit OK",1.0);
                h_fitprobability.Fill(fit_result.Probability);
                utils::ParticleTools::FillIMCombinations([this] (double x) {IM_2g_byFit.Fill(x);},  2, photons);
                utils::ParticleTools::FillIMCombinations([this] (double x) {IM_2g_fitted.Fill(x);},  2, fitter.GetFittedPhotons());
            }
            else {
                steps->Fill("Fit failed",1.0);
            }

        }
    }
}

void KinFitPi0::MultiPi0::ShowResult()
{
    canvas(std_ext::formatter() << "KinFitPi0: " << multiplicity << "Pi0")
            << steps
            << Proton_Coplanarity
            << h_taggtime
            << h_missingmass.subtracted
            << h_fitprobability.subtracted
            << IM_2g_byFit.subtracted
            << IM_2g_fitted.subtracted
            << endc;
}

AUTO_REGISTER_PHYSICS(KinFitPi0)
