#include "checkkinfit.h"

#include "expconfig/ExpConfig.h"
#include "utils/uncertainties/Interpolated.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "utils/ProtonPhotonCombs.h"
#include "utils/ParticleTools.h"
#include "utils/Combinatorics.h"

#include "base/Logger.h"

#include "TLorentzVector.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

scratch_lheijken_checkkinfit::scratch_lheijken_checkkinfit(const std::string& name, OptionsPtr opts) :
    Physics(name, opts),
    promptrandom(ExpConfig::Setup::Get()),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                utils::UncertaintyModels::Interpolated::Type_t::MC,
                make_shared<utils::UncertaintyModels::FitterSergey>())),
    fitter(nullptr, opts->Get<bool>("FitZVertex", true),MakeFitSettings(30)),
    npfitter(nullptr, opts->Get<bool>("FitZVertex", true),MakeFitSettings(30))
{
    fitter.SetZVertexSigma(3.0);
    npfitter.SetZVertexSigma(3.0);

    CreateHistos();

    taps_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
    TAPSZPos = taps_detector->GetZPosition();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    CBRad = cb_detector->GetInnerRadius();
}

void scratch_lheijken_checkkinfit::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    if(!triggersimu.HasTriggered())
        return;
    h_Steps->Fill(0);

    fitter.SetUncertaintyModel(fit_model);
    npfitter.SetUncertaintyModel(fit_model);

    const TEventData& data = event.Reconstructed();

    //-- Check if the data is MC and classify which decay channel is analysed
    //   When analysing MC events, both geant and pluto input are required
    bool isMC=false;
    string decay;
    if(data.ID.isSet(ant::TID::Flags_t::MC) && event.MCTrue().ParticleTree){
        decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
        isMC=true;
    }
    else
        decay = "data" + ExpConfig::Setup::Get().GetName();
    vector<bool> WhichMC; // 0 : pi0->eeg, 1 : pi0->gg, 2 : pi0pi0->gggg
    bool MCpi0eeg = false; bool MCpi0gg = false; bool MC2pi04g = false;
    double true_z_vert = 0;
    if(isMC){
        MCpi0eeg  = event.MCTrue().ParticleTree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_eeg),
                                                         utils::ParticleTools::MatchByParticleName);
        MCpi0gg  = event.MCTrue().ParticleTree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),
                                                        utils::ParticleTools::MatchByParticleName);
        MC2pi04g  = event.MCTrue().ParticleTree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::TwoPi0_4g),
                                                        utils::ParticleTools::MatchByParticleName);
        const TEventData& mctrue = event.MCTrue();
        true_z_vert = mctrue.Target.Vertex.z;
    }
    WhichMC.push_back(MCpi0eeg); WhichMC.push_back(MCpi0gg); WhichMC.push_back(MC2pi04g);

    //-- The reconstructed candidates
    const auto& candidates = data.Candidates;
    int nrcands = candidates.size();
    //-- At the moment, atleast the number of final state particles, exluding the proton
    //    pi0->eeg : 3,  pi0->gg : 2, pi0pi0->gggg : 4
    int nrphotonsinfit = 3;
    //-- Also at the moment, when analysing MC, let the decay type decide
    if(MCpi0eeg){ nrphotonsinfit = 3; }
    if(MCpi0gg){ nrphotonsinfit = 2; }
    if(MC2pi04g){ nrphotonsinfit = 4; }
    if(nrcands < nrphotonsinfit) return;
    h_Steps->Fill(1);

    //-- Make all possible combinations of the candidates into groups containing one proton and the rest photons
    utils::ProtonPhotonCombs proton_photons(candidates);
    //-- Also make all possible combinations of nrphotonsinfit

    //-- MC stuff, for the pi0 -> eeg,  pi0->gg and pi0pi0->4g channels
    vector<TParticlePtr> TruePart;
    vector<vector<TParticlePtr>> TrueRecoMatchPart;
    if(isMC){
        //-- Fetch the true final state particles
        if(MCpi0eeg){
            auto TruePartp_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton,event.MCTrue().ParticleTree);
            auto TruePartep_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::ePlus,event.MCTrue().ParticleTree);
            auto TruePartem_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::eMinus,event.MCTrue().ParticleTree);
            auto TruePartg_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Photon,event.MCTrue().ParticleTree);
            TruePart.push_back(TruePartp_eeg); TruePart.push_back(TruePartep_eeg);
            TruePart.push_back(TruePartem_eeg); TruePart.push_back(TruePartg_eeg);
            h_IMeeg_True->Fill((*TruePartep_eeg + *TruePartem_eeg + *TruePartg_eeg).M());
        }
        if(MCpi0gg){
            auto TruePartp_gg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton,event.MCTrue().ParticleTree);
            auto TruePartgs_gg = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, event.MCTrue().ParticleTree);
            if(TruePartgs_gg.size() != 2){LOG(INFO)<<"Didn't find two photons in pi0->gg trueMC"; return;}
            TruePart.push_back(TruePartp_gg);
            TruePart.push_back(TruePartgs_gg.front()); TruePart.push_back(TruePartgs_gg.back());
            h_IMgg_True->Fill((*TruePartgs_gg.front() + *TruePartgs_gg.back()).M());
        }
        if(MC2pi04g){
            auto TruePartp_gg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton,event.MCTrue().ParticleTree);
            auto TruePartgs_gg = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, event.MCTrue().ParticleTree);
            if(TruePartgs_gg.size() != 4){LOG(INFO)<<"Didn't find four photons in pi0pi0->gggg trueMC"; return;}
            TruePart.push_back(TruePartp_gg);
            TruePart.push_back(TruePartgs_gg.at(0)); TruePart.push_back(TruePartgs_gg.at(1));
            TruePart.push_back(TruePartgs_gg.at(2)); TruePart.push_back(TruePartgs_gg.at(3));
            h_IMgg_True->Fill((*TruePartgs_gg.at(0) + *TruePartgs_gg.at(1)).M());
            h_IMgg_True->Fill((*TruePartgs_gg.at(2) + *TruePartgs_gg.at(3)).M());
        }
        //-- Fetch the true beamtarget particle
        auto TruePart_beam = utils::ParticleTools::FindParticle(ParticleTypeDatabase::BeamTarget,event.MCTrue().ParticleTree);
        TruePart.push_back(TruePart_beam);

        //-- Match the true particles to the reconstructed candidates
        auto mcparticleslist = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
        const auto mcparticles = mcparticleslist.GetAll();
        DoMatchTrueRecoStuff(mcparticles, TruePart, candidates, TrueRecoMatchPart, true_z_vert);

        //---TEMP---
        if(chooserecprot>0){
            bool protonrecmatch = false;
            for(int j=0; j<(int)TrueRecoMatchPart.size(); j++){
                TParticlePtr tmph = TrueRecoMatchPart.at(j).at(0);
                if(tmph->Type() == ParticleTypeDatabase::Proton) protonrecmatch = true;
            }
            if(chooserecprot==1 && !protonrecmatch) return;
            if(chooserecprot==2 && protonrecmatch) return;
        }
    }

    //-- Tagger loop
    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        //-- Only bother with tagger hits which are inside the prompt or random windows
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        double taggw = promptrandom.FillWeight();
        h_Steps->Fill(2,taggw);
        TLorentzVector InitialPhotonVec = taggerhit.GetPhotonBeam();

        //-- compare Ebeam True-Rec
        TParticlePtr truebeam = 0;
        for(auto& truepart: TruePart){
            if(truepart->Type() == ParticleTypeDatabase::BeamTarget){
                truebeam = truepart;
                h_TrueRec_Ebeam->Fill(taggerhit.PhotonEnergy,truebeam->Ek()-taggerhit.PhotonEnergy,taggw);
                break;
            }
        }


        TParticleList recphotonsbestfit_measp; vector<LorentzVec> fitphotonsbestfit_measp;
        TParticlePtr recprotonbestfit_measp = 0; LorentzVec fitprotonbestfit_measp;
        std::vector<utils::Fitter::FitParticle> fitparticlesbestfit_measp;
        double fitprob_measp = -50;
        double fit_z_vert_measp = -50;
        double fitbeamE_measp = -50;
        double pullbeamE_measp = -50;
        double pullzvert_measp = -50;

        //-- Do a kinematic fit for 1 proton (only energy unmeasured) and a fixed number of photons -- currently only if at least #cand >= # all f.s. particles
        if(nrcands > nrphotonsinfit){
            h_Steps->Fill(3,taggw);
            //-- create the 1 proton - X photons combinations, require the MM(Xg) to be close to proton mass
            auto filtered_combs = proton_photons()
                                  .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(350).Round());
            if(filtered_combs.empty()) {
                continue;
            }
            //-- loop over the (filtered) proton combinations
            for(const auto& comb : filtered_combs) {
                TParticlePtr protontofit = comb.Proton;
                //-- select on the number of photons required in the current kinfit hypothesis
                for(auto photoncomb = analysis::utils::makeCombination(comb.Photons,nrphotonsinfit); !photoncomb.done(); ++photoncomb ){
                    TParticleList photonstofit;
                    for(int i=0; i<nrphotonsinfit; i++)
                        photonstofit.push_back(photoncomb.at(i));

                    const auto& result = fitter.DoFit(taggerhit.PhotonEnergy, protontofit, photonstofit);

                    //-- If the fit did not converge or did not result in a better probability, try next combination
                    if(result.Status != APLCON::Result_Status_t::Success)
                        continue;
                    if(fitprob_measp > result.Probability)
                        continue;
                    //-- Else, store the fit result and the reconstructed input
                    fitprob_measp = result.Probability;
                    recphotonsbestfit_measp = photonstofit;
                    fitphotonsbestfit_measp.clear();
                    for(const auto& fph : fitter.GetFittedPhotons())
                        fitphotonsbestfit_measp.push_back(*fph);
                    recprotonbestfit_measp = protontofit;
                    fitprotonbestfit_measp = *fitter.GetFittedProton();
                    fit_z_vert_measp = fitter.GetFittedZVertex();
                    fitbeamE_measp = fitter.GetFittedBeamE();
                    fitparticlesbestfit_measp.clear();
                    fitparticlesbestfit_measp = fitter.GetFitParticles();
                    pullbeamE_measp = fitter.GetBeamEPull();
                    pullzvert_measp = fitter.GetZVertexPull();
                }
            }
            //---TEMP----
            bool keep = true;
            if(chooseimreg>0){
                double imggtrh = 0;
                if(fitphotonsbestfit_measp.size()>0){
                    LorentzVec phsum;
                    for(const auto& ph: fitphotonsbestfit_measp)
                        phsum += ph;
                    imggtrh = phsum.M();
                }
                if(chooseimreg==1 && (imggtrh < 100 || imggtrh > 170)) keep = false;
                if(chooseimreg==2 && imggtrh<200) keep = false;
            }
            if(keep){
            //-- Make cuts on the fit probability
            int kfnr = 0;
            double pcut[] = {0.,0.001,0.01};
            for(int ipc=0; ipc<nrPcut; ipc++){
                if(fitprob_measp >= pcut[ipc]) {
                    h_Probability[kfnr][ipc]->Fill(fitprob_measp,taggw);
                    //-- Make all the comparisons between the protons and photons
                    DoFitComparisons(kfnr,ipc,nrphotonsinfit,fitprotonbestfit_measp,recprotonbestfit_measp,fitphotonsbestfit_measp,recphotonsbestfit_measp,TrueRecoMatchPart,TruePart,fitparticlesbestfit_measp,fit_z_vert_measp,taggw);
                    //-- Also for the beam and z-vertex
                    if(truebeam){
                        h_TrueFit_Ebeam[kfnr][ipc]->Fill(truebeam->Ek(),truebeam->Ek() - fitbeamE_measp,taggw);
                        h_TrueFit_zvert[kfnr][ipc]->Fill(true_z_vert,true_z_vert - fit_z_vert_measp,taggw);
                    }
                    h_FitRec_Ebeam[kfnr][ipc]->Fill(taggerhit.PhotonEnergy, fitbeamE_measp - taggerhit.PhotonEnergy,taggw);
                    h_EbeamPulls[kfnr][ipc]->Fill(pullbeamE_measp,taggw);
                    h_ZvertPulls[kfnr][ipc]->Fill(pullzvert_measp,taggw);
                }
            }
            }
        }

        TParticleList recphotonsbestfit_unmeasp; vector<LorentzVec> fitphotonsbestfit_unmeasp;
        TParticlePtr recprotonbestfit_unmeasp = 0; LorentzVec fitprotonbestfit_unmeasp;
        std::vector<utils::Fitter::FitParticle> fitparticlesbestfit_unmeasp;
        double fitprob_unmeasp = -50;
        double fit_z_vert_unmeasp = -50;
        double fitbeamE_unmeasp = -50;
        double pullbeamE_unmeasp = -50;
        double pullzvert_unmeasp = -50;

        //-- Do a fit with the proton unmeasured and a fixed number of photons  -- currently only if at least #cand >= # the required number of photons
        if(nrcands >= nrphotonsinfit){
            h_Steps->Fill(4,taggw);
            //-- Make combinations of all photon candidates -- if their missing mass is near the proton mass
            TParticleList photonlist;
            for(auto cand : candidates.get_iter())
                photonlist.emplace_back(std::make_shared<TParticle>(ParticleTypeDatabase::Photon,cand));
            for(auto photoncomb = analysis::utils::makeCombination(photonlist,nrphotonsinfit); !photoncomb.done(); ++photoncomb ){
                //-- check that the MM(p) is near mp
                TLorentzVector photonssum;
                for(int ipc=0; ipc<nrphotonsinfit; ipc++)
                    photonssum += *photoncomb.at(ipc);
                double mmp = (LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass()) + InitialPhotonVec - photonssum).M();
                if(!ParticleTypeDatabase::Proton.GetWindow(350).Round().Contains(mmp))
                    continue;
                TParticleList photonstofit;
                for(int i=0; i<nrphotonsinfit; i++)
                    photonstofit.push_back(photoncomb.at(i));

                const auto& result = npfitter.DoFit(taggerhit.PhotonEnergy, photonstofit);

                //-- If the fit did not converge or did not result in a better probability, try next combination
                if(result.Status != APLCON::Result_Status_t::Success)
                    continue;
                if(fitprob_unmeasp > result.Probability)
                    continue;
                //-- Else, store the fit result and the reconstructed input
                fitprob_unmeasp = result.Probability;
                recphotonsbestfit_unmeasp = photonstofit;
                fitphotonsbestfit_unmeasp.clear();
                for(const auto& fph : npfitter.GetFittedPhotons())
                    fitphotonsbestfit_unmeasp.push_back(*fph);
                recprotonbestfit_unmeasp = 0;
                fitprotonbestfit_unmeasp = npfitter.GetFittedProton();
                fit_z_vert_unmeasp = npfitter.GetFittedZVertex();
                fitbeamE_unmeasp = npfitter.GetFittedBeamE();
                fitparticlesbestfit_unmeasp.clear();
                fitparticlesbestfit_unmeasp = npfitter.GetFitParticles();
                pullbeamE_unmeasp = npfitter.GetBeamEPull();
                pullzvert_unmeasp = npfitter.GetZVertexPull();
            }

            //---TEMP----
            bool keep = true;
            if(chooseimreg>0){
                double imggtrh = 0;
                if(fitphotonsbestfit_unmeasp.size()>0){
                    LorentzVec phsum;
                    for(const auto& ph: fitphotonsbestfit_unmeasp)
                        phsum += ph;
                    imggtrh = phsum.M();
                }
                if(chooseimreg==1 && (imggtrh < 100 || imggtrh > 170)) keep = false;
                if(chooseimreg==2 && imggtrh<200) keep = false;
            }
            if(keep){
            //-- Make cuts on the fit probability
            int kfnr = 1;
            double pcut[] = {0.,0.001,0.01};
            for(int ipc=0; ipc<nrPcut; ipc++){
                if(fitprob_unmeasp >= pcut[ipc]) {
                    h_Probability[kfnr][ipc]->Fill(fitprob_unmeasp,taggw);
                    //-- Make all the comparisons between the protons and photons
                    DoFitComparisons(kfnr,ipc,nrphotonsinfit,fitprotonbestfit_unmeasp,recprotonbestfit_unmeasp,fitphotonsbestfit_unmeasp,recphotonsbestfit_unmeasp,TrueRecoMatchPart,TruePart,fitparticlesbestfit_unmeasp,fit_z_vert_unmeasp,taggw);
                    //-- Also for the beam and z-vertex
                    if(truebeam){
                        h_TrueFit_Ebeam[kfnr][ipc]->Fill(truebeam->Ek(),truebeam->Ek() - fitbeamE_unmeasp,taggw);
                        h_TrueFit_zvert[kfnr][ipc]->Fill(true_z_vert,true_z_vert - fit_z_vert_unmeasp,taggw);
                    }
                    h_FitRec_Ebeam[kfnr][ipc]->Fill(taggerhit.PhotonEnergy, fitbeamE_unmeasp - taggerhit.PhotonEnergy,taggw);
                    h_EbeamPulls[kfnr][ipc]->Fill(pullbeamE_unmeasp,taggw);
                    h_ZvertPulls[kfnr][ipc]->Fill(pullzvert_unmeasp,taggw);
                }
            }
            }
        }

        //-- Compare the two fit results and store histograms in two additional cases
        if(fitprob_measp >=0) h_Steps->Fill(5,taggw);
        if(fitprob_unmeasp >=0) h_Steps->Fill(6,taggw);
        //-- Measured proton fit is succesful and more (or equally) probable
        if((fitprob_measp >= 0) && (fitprob_measp >= fitprob_unmeasp)){
            if((fitbeamE_measp == fitbeamE_unmeasp)) h_Steps->Fill(9,taggw);
            else h_Steps->Fill(7,taggw);
            //---TEMP----
            bool keep = true;
            if(chooseimreg>0){
                double imggtrh = 0;
                if(fitphotonsbestfit_measp.size()>0){
                    LorentzVec phsum;
                    for(const auto& ph: fitphotonsbestfit_measp)
                        phsum += ph;
                    imggtrh = phsum.M();
                }
                if(chooseimreg==1 && (imggtrh < 100 || imggtrh > 170)) keep = false;
                if(chooseimreg==2 && imggtrh<200) keep = false;
            }
            if(keep){
            //-- Make cuts on the fit probability
            int kfnr = 2; int kfnrboth = 4;
            double pcut[] = {0.,0.001,0.01};
            for(int ipc=0; ipc<nrPcut; ipc++){
                if(fitprob_measp >= pcut[ipc]) {
                    h_Probability[kfnr][ipc]->Fill(fitprob_measp,taggw);
                    h_Probability[kfnrboth][ipc]->Fill(fitprob_measp,taggw);
                    //-- Make all the comparisons between the protons and photons
                    DoFitComparisons(kfnr,ipc,nrphotonsinfit,fitprotonbestfit_measp,recprotonbestfit_measp,fitphotonsbestfit_measp,recphotonsbestfit_measp,TrueRecoMatchPart,TruePart,fitparticlesbestfit_measp,fit_z_vert_measp,taggw);
                    DoFitComparisons(kfnrboth,ipc,nrphotonsinfit,fitprotonbestfit_measp,recprotonbestfit_measp,fitphotonsbestfit_measp,recphotonsbestfit_measp,TrueRecoMatchPart,TruePart,fitparticlesbestfit_measp,fit_z_vert_measp,taggw);
                    //-- Also for the beam and z-vertex
                    if(truebeam){
                        h_TrueFit_Ebeam[kfnr][ipc]->Fill(truebeam->Ek(),truebeam->Ek() - fitbeamE_measp,taggw);
                        h_TrueFit_zvert[kfnr][ipc]->Fill(true_z_vert,true_z_vert - fit_z_vert_measp,taggw);
                        h_TrueFit_Ebeam[kfnrboth][ipc]->Fill(truebeam->Ek(),truebeam->Ek() - fitbeamE_measp,taggw);
                        h_TrueFit_zvert[kfnrboth][ipc]->Fill(true_z_vert,true_z_vert - fit_z_vert_measp,taggw);
                    }
                    h_FitRec_Ebeam[kfnr][ipc]->Fill(taggerhit.PhotonEnergy, fitbeamE_measp - taggerhit.PhotonEnergy,taggw);
                    h_FitRec_Ebeam[kfnrboth][ipc]->Fill(taggerhit.PhotonEnergy, fitbeamE_measp - taggerhit.PhotonEnergy,taggw);
                    h_EbeamPulls[kfnr][ipc]->Fill(pullbeamE_measp,taggw);
                    h_EbeamPulls[kfnrboth][ipc]->Fill(pullbeamE_measp,taggw);
                    h_ZvertPulls[kfnr][ipc]->Fill(pullzvert_measp,taggw);
                    h_ZvertPulls[kfnrboth][ipc]->Fill(pullzvert_measp,taggw);
                }
            }
            }
        }
        //-- Unmeasured proton fit is succsessful and more probable
        if((fitprob_unmeasp >= 0) && (fitprob_unmeasp > fitprob_measp)){
            h_Steps->Fill(8,taggw);
            //---TEMP----
            bool keep = true;
            if(chooseimreg>0){
                double imggtrh = 0;
                if(fitphotonsbestfit_unmeasp.size()>0){
                    LorentzVec phsum;
                    for(const auto& ph: fitphotonsbestfit_unmeasp)
                        phsum += ph;
                    imggtrh = phsum.M();
                }
                if(chooseimreg==1 && (imggtrh < 100 || imggtrh > 170)) keep = false;
                if(chooseimreg==2 && imggtrh<200) keep = false;
            }
            if(keep){
            //-- Make cuts on the fit probability
            int kfnr = 3; int kfnrboth = 4;
            double pcut[] = {0.,0.001,0.01};
            for(int ipc=0; ipc<nrPcut; ipc++){
                if(fitprob_unmeasp >= pcut[ipc]) {
                    h_Probability[kfnr][ipc]->Fill(fitprob_unmeasp,taggw);
                    h_Probability[kfnrboth][ipc]->Fill(fitprob_unmeasp,taggw);
                    //-- Make all the comparisons between the protons and photons
                    DoFitComparisons(kfnr,ipc,nrphotonsinfit,fitprotonbestfit_unmeasp,recprotonbestfit_unmeasp,fitphotonsbestfit_unmeasp,recphotonsbestfit_unmeasp,TrueRecoMatchPart,TruePart,fitparticlesbestfit_unmeasp,fit_z_vert_unmeasp,taggw);
                    DoFitComparisons(kfnrboth,ipc,nrphotonsinfit,fitprotonbestfit_unmeasp,recprotonbestfit_unmeasp,fitphotonsbestfit_unmeasp,recphotonsbestfit_unmeasp,TrueRecoMatchPart,TruePart,fitparticlesbestfit_unmeasp,fit_z_vert_unmeasp,taggw);
                    //-- Also for the beam and z-vertex
                    if(truebeam){
                        h_TrueFit_Ebeam[kfnr][ipc]->Fill(truebeam->Ek(),truebeam->Ek() - fitbeamE_unmeasp,taggw);
                        h_TrueFit_zvert[kfnr][ipc]->Fill(true_z_vert,true_z_vert - fit_z_vert_unmeasp,taggw);
                        h_TrueFit_Ebeam[kfnrboth][ipc]->Fill(truebeam->Ek(),truebeam->Ek() - fitbeamE_unmeasp,taggw);
                        h_TrueFit_zvert[kfnrboth][ipc]->Fill(true_z_vert,true_z_vert - fit_z_vert_unmeasp,taggw);
                    }
                    h_FitRec_Ebeam[kfnr][ipc]->Fill(taggerhit.PhotonEnergy, fitbeamE_unmeasp - taggerhit.PhotonEnergy,taggw);
                    h_FitRec_Ebeam[kfnrboth][ipc]->Fill(taggerhit.PhotonEnergy, fitbeamE_unmeasp - taggerhit.PhotonEnergy,taggw);
                    h_EbeamPulls[kfnr][ipc]->Fill(pullbeamE_unmeasp,taggw);
                    h_EbeamPulls[kfnrboth][ipc]->Fill(pullbeamE_unmeasp,taggw);
                    h_ZvertPulls[kfnr][ipc]->Fill(pullzvert_unmeasp,taggw);
                    h_ZvertPulls[kfnrboth][ipc]->Fill(pullzvert_unmeasp,taggw);
                }
            }
            }
        }
    }
}

void scratch_lheijken_checkkinfit::ShowResult()
{
    canvas(GetName())
            << h_Steps
            << endc;
}

void scratch_lheijken_checkkinfit::CreateHistos()
{
    //-- Histogram folders
    auto hfOverview = new HistogramFactory("Overview",HistFac,"");
    auto hfTrueRec = new HistogramFactory("TrueRec",HistFac,"");
    auto hfTrueRecCB = new HistogramFactory("CB",*hfTrueRec);
    auto hfTrueRecTA = new HistogramFactory("TAPS",*hfTrueRec);
    auto hfFitRec = new HistogramFactory("FitRec",HistFac,"");
    auto hfTrueFit = new HistogramFactory("TrueFit",HistFac,"");
    auto hfPulls = new HistogramFactory("Pulls",HistFac,"");

    //-- Histogram binnings and names
    string varname[] = {"Ek","theta","phi"};
    string vartitle[] = {"Ek","#theta","#phi"};
    string partname[] = {"p","ep","em","g"};
    string parttitle[] = {"p","e^{+}","e^{-}","#gamma"};
    string kfname[] = {"measp_all","unmeasp_all","measp_sel","unmeasp_sel","both_sel"};
    string kftitle[] = {"kf measured p all", "kf unmeasured p all","kf measured p sel", "kf unmeasured p sel","kf both sel"};
    string pcutname[] = {"pc0","pc001","pc01"};
    string pcuttitle[] = {"P(#chi^{2})>0","P(#chi^{2})>0.001","P(#chi^{2})>0.01"};
    string fitvarnameCB[] = {"invEk","theta","phi","R"};
    string fitvarnameTA[] = {"invEk","Rxy","phi","L"};
    string fitvartitleCB[] = {"1/Ek","#theta","#phi","R"};
    string fitvartitleTA[] = {"1/Ek","R_{xy}","#phi","L"};
    const BinSettings EkBins = BinSettings(500,0,1000);
    const BinSettings ThCBBins = BinSettings(90,0,180);
    const BinSettings ThTABins = BinSettings(40,0,40);
    const BinSettings PhBins = BinSettings(200,-200,200);
    vector<BinSettings> VarBins;
    VarBins.push_back(BinSettings(101,-100,100));
    VarBins.push_back(BinSettings(101,-100,100));
    VarBins.push_back(BinSettings(101,-100,100));

    //-- The histograms
    h_Steps = hfOverview->makeTH1D("Steps","","",BinSettings(15),"h_Steps");
    h_Steps->GetXaxis()->SetBinLabel(1,"CBEsum"); h_Steps->GetXaxis()->SetBinLabel(2,">=nrphot"); h_Steps->GetXaxis()->SetBinLabel(3,"PR tagghit");
    h_Steps->GetXaxis()->SetBinLabel(4,">nrphot"); h_Steps->GetXaxis()->SetBinLabel(5,">=nrphot");
    h_Steps->GetXaxis()->SetBinLabel(6,"measpSuc"); h_Steps->GetXaxis()->SetBinLabel(7,"unmeaspSuc");
    h_Steps->GetXaxis()->SetBinLabel(8,"measpPlar"); h_Steps->GetXaxis()->SetBinLabel(9,"unmeaspPlarg"); h_Steps->GetXaxis()->SetBinLabel(10,"Peq");
    for(int i=0; i<nrFitCases; i++){
        //-- Create subfolders for each case
        auto hfOverviewKFCasestemp = new HistogramFactory(kfname[i],*hfOverview);
        auto hfFitRecKFCasestemp = new HistogramFactory(kfname[i],*hfFitRec);
        auto hfTrueFitKFCasestemp = new HistogramFactory(kfname[i],*hfTrueFit);
        auto hfPullsKFCasestemp = new HistogramFactory(kfname[i],*hfPulls);
        for(int j=0; j<nrPcut; j++){
            //-- Create subfolders for each cut
            auto hfOverviewPcutstemp = new HistogramFactory(pcutname[j],*hfOverviewKFCasestemp);
            auto hfFitRecPcutstemp = new HistogramFactory(pcutname[j],*hfFitRecKFCasestemp);
            auto hfFitRecPcutstempCB = new HistogramFactory("CB",*hfFitRecPcutstemp);
            auto hfFitRecPcutstempTA = new HistogramFactory("TAPS",*hfFitRecPcutstemp);
            auto hfTrueFitPcutstemp = new HistogramFactory(pcutname[j], *hfTrueFitKFCasestemp);
            auto hfTrueFitPcutstempCB = new HistogramFactory("CB", *hfTrueFitPcutstemp);
            auto hfTrueFitPcutstempTA = new HistogramFactory("TAPS", *hfTrueFitPcutstemp);
            auto hfPullsPcutstemp = new HistogramFactory(pcutname[j], *hfPullsKFCasestemp);
            auto hfPullsPcutstempCB = new HistogramFactory("CB", *hfPullsPcutstemp);
            auto hfPullsPcutstempTA = new HistogramFactory("TAPS", *hfPullsPcutstemp);

            h_Probability[i][j] = hfOverviewPcutstemp->makeTH1D(Form("Fit probability, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"P(#chi^{2})","",BinSettings(1000,0.,1.),Form("h_Probability_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            //--- The kinematical variables for the final state particles
            for(int k=0; k<nrPartTypes; k++){
                for(int l=0; l<nrKinVars; l++){
                    if(i==0 && j==0){
                        //-- TrueRec
                        h_TrueRec_vsEk_CB[k][l] = hfTrueRecCB->makeTH2D(Form("TrueRec %s vs Ek^{rec} CB %s",vartitle[l].c_str(),parttitle[k].c_str()),"Ek [MeV]",Form("True - Rec %s",vartitle[l].c_str()),EkBins,VarBins.at(l),Form("hTrueRec%svsEk_CB_%s",varname[l].c_str(),partname[k].c_str()),true);
                        h_TrueRec_vsEk_TAPS[k][l] = hfTrueRecTA->makeTH2D(Form("TrueRec %s vs Ek^{rec} TAPS %s",vartitle[l].c_str(),parttitle[k].c_str()),"Ek [MeV]",Form("True - Rec %s",vartitle[l].c_str()),EkBins,VarBins.at(l),Form("hTrueRec%svsEk_TAPS_%s",varname[l].c_str(),partname[k].c_str()),true);
                        h_TrueRec_vsTh_CB[k][l] = hfTrueRecCB->makeTH2D(Form("TrueRec %s vs Th^{rec} CB %s",vartitle[l].c_str(),parttitle[k].c_str()),"#theta [deg]",Form("True - Rec %s",vartitle[l].c_str()),ThCBBins,VarBins.at(l),Form("hTrueRec%svsTh_CB_%s",varname[l].c_str(),partname[k].c_str()),true);
                        h_TrueRec_vsTh_TAPS[k][l] = hfTrueRecTA->makeTH2D(Form("TrueRec %s vs Th^{rec} TAPS %s",vartitle[l].c_str(),parttitle[k].c_str()),"#theta [deg]",Form("True - Rec %s",vartitle[l].c_str()),ThTABins,VarBins.at(l),Form("hTrueRec%svsTh_TAPS_%s",varname[l].c_str(),partname[k].c_str()),true);
                        h_TrueRec_vsPh_CB[k][l] = hfTrueRecCB->makeTH2D(Form("TrueRec %s vs Ph^{rec} CB %s",vartitle[l].c_str(),parttitle[k].c_str()),"#phi [deg]",Form("True - Rec %s",vartitle[l].c_str()),PhBins,VarBins.at(l),Form("hTrueRec%svsPh_CB_%s",varname[l].c_str(),partname[k].c_str()),true);
                        h_TrueRec_vsPh_TAPS[k][l] = hfTrueRecTA->makeTH2D(Form("TrueRec %s vs Ph^{rec} TAPS %s",vartitle[l].c_str(),parttitle[k].c_str()),"#phi [deg]",Form("True - Rec %s",vartitle[l].c_str()),PhBins,VarBins.at(l),Form("hTrueRec%svsPh_TAPS_%s",varname[l].c_str(),partname[k].c_str()),true);
                    }
                    //-- FitRec
                    h_FitRec_vsEk_CB[i][j][k][l] = hfFitRecPcutstempCB->makeTH2D(Form("FitRec %s vs Ek^{rec} CB %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"Ek [MeV]",Form("Fit - Rec %s",vartitle[l].c_str()),EkBins,VarBins.at(l),Form("hFitRec%svsEk_CB_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_FitRec_vsEk_TAPS[i][j][k][l] = hfFitRecPcutstempTA->makeTH2D(Form("FitRec %s vs Ek^{rec} TAPS %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"Ek [MeV]",Form("Fit - Rec %s",vartitle[l].c_str()),EkBins,VarBins.at(l),Form("hFitRec%svsEk_TAPS_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_FitRec_vsTh_CB[i][j][k][l] = hfFitRecPcutstempCB->makeTH2D(Form("FitRec %s vs Th^{rec} CB %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"#theta [deg]",Form("Fit - Rec %s",vartitle[l].c_str()),ThCBBins,VarBins.at(l),Form("hFitRec%svsTh_CB_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_FitRec_vsTh_TAPS[i][j][k][l] = hfFitRecPcutstempTA->makeTH2D(Form("FitRec %s vs Th^{rec} TAPS %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"#theta [deg]",Form("Fit - Rec %s",vartitle[l].c_str()),ThTABins,VarBins.at(l),Form("hFitRec%svsTh_TAPS_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_FitRec_vsPh_CB[i][j][k][l] = hfFitRecPcutstempCB->makeTH2D(Form("FitRec %s vs Ph^{rec} CB %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"#phi [deg]",Form("Fit - Rec %s",vartitle[l].c_str()),PhBins,VarBins.at(l),Form("hFitRec%svsPh_CB_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_FitRec_vsPh_TAPS[i][j][k][l] = hfFitRecPcutstempTA->makeTH2D(Form("FitRec %s vs Ph^{rec} TAPS %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"#phi [deg]",Form("Fit - Rec %s",vartitle[l].c_str()),PhBins,VarBins.at(l),Form("hFitRec%svsPh_TAPS_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    //-- TrueFit
                    h_TrueFit_vsEk_CB[i][j][k][l] = hfTrueFitPcutstempCB->makeTH2D(Form("TrueFit %s vs Ek^{true} CB %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"Ek [MeV]",Form("True - Fit %s",vartitle[l].c_str()),EkBins,VarBins.at(l),Form("hTrueFit%svsEk_CB_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_TrueFit_vsEk_TAPS[i][j][k][l] = hfTrueFitPcutstempTA->makeTH2D(Form("TrueFit %s vs Ek^{true} TAPS %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"Ek [MeV]",Form("True - Fit %s",vartitle[l].c_str()),EkBins,VarBins.at(l),Form("hTrueFit%svsEk_TAPS_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_TrueFit_vsTh_CB[i][j][k][l] = hfTrueFitPcutstempCB->makeTH2D(Form("TrueFit %s vs Th^{true} CB %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"#theta [deg]",Form("True - Fit %s",vartitle[l].c_str()),ThCBBins,VarBins.at(l),Form("hTrueFit%svsTh_CB_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_TrueFit_vsTh_TAPS[i][j][k][l] = hfTrueFitPcutstempTA->makeTH2D(Form("TrueFit %s vs Th^{true} TAPS %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"#theta [deg]",Form("True - Fit %s",vartitle[l].c_str()),ThTABins,VarBins.at(l),Form("hTrueFit%svsTh_TAPS_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_TrueFit_vsPh_CB[i][j][k][l] = hfTrueFitPcutstempCB->makeTH2D(Form("TrueFit %s vs Ph^{true} CB %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"#phi [deg]",Form("True - Fit %s",vartitle[l].c_str()),PhBins,VarBins.at(l),Form("hTrueFit%svsPh_CB_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_TrueFit_vsPh_TAPS[i][j][k][l] = hfTrueFitPcutstempTA->makeTH2D(Form("TrueFit %s vs Ph^{true} TAPS %s, %s, %s",vartitle[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),"#phi [deg]",Form("True - Fit %s",vartitle[l].c_str()),PhBins,VarBins.at(l),Form("hTrueFit%svsPh_TAPS_%s_%s_%s",varname[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                }
                for(int l=0; l<nrFitVars; l++){
                    //-- Pulls
                    h_PartPulls_CB[i][j][k][l] = hfPullsPcutstempCB->makeTH1D(Form("%s pull for %s, %s, %s",fitvartitleCB[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),Form("%s pull",fitvartitleCB[l].c_str()),"",BinSettings(200,-10.,10.),Form("h_PartPulls_%s_%s_%s_%s",fitvarnameCB[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                    h_PartPulls_TAPS[i][j][k][l] = hfPullsPcutstempTA->makeTH1D(Form("%s pull for %s, %s, %s",fitvartitleTA[l].c_str(),parttitle[k].c_str(),kftitle[i].c_str(),pcuttitle[j].c_str()),Form("%s pull",fitvartitleTA[l].c_str()),"",BinSettings(200,-10.,10.),Form("h_PartPulls_%s_%s_%s_%s",fitvarnameTA[l].c_str(),partname[k].c_str(),kfname[i].c_str(),pcutname[j].c_str()),true);
                }
            }
            //--- The Ebeam and z-vertex histos
            h_TrueFit_zvert[i][j] = hfTrueFitPcutstemp->makeTH2D(Form("True-Fit z-vertex vs True, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"z_{true}","z_{true} - z_{fit}",BinSettings(50,-15,15),BinSettings(51,-10,10),Form("h_TrueFit_zvert_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            if(i==0 && j==0) h_TrueRec_Ebeam = hfTrueRec->makeTH2D("TrueRec E_{beam} vs E_{beam}^{rec}","E_{beam}^{rec}","True - Rec E_{beam}",BinSettings(1600,0.,1600.),BinSettings(101,-10,10),"h_TrueRec_Ebeam",true);
            h_FitRec_Ebeam[i][j] = hfFitRecPcutstemp->makeTH2D(Form("FitRec E_{beam} vs E_{beam}^{rec}, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"E_{beam}^{rec}","Fit - Rec E_{beam}",BinSettings(1600,0.,1600.),BinSettings(101,-10,10),Form("h_FitRec_Ebeam_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            h_TrueFit_Ebeam[i][j] = hfTrueFitPcutstemp->makeTH2D(Form("TrueFit E_{beam} vs E_{beam}^{true}, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"E_{beam}^{true}","True - Fit E_{beam}",BinSettings(1600,0.,1600.),BinSettings(101,-10,10),Form("h_TrueFit_Ebeam_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            h_EbeamPulls[i][j] = hfPullsPcutstemp->makeTH1D(Form("Ebeam pull, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"Pull","",BinSettings(200,-10.,10.),Form("h_EbeamPulls_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            h_ZvertPulls[i][j] = hfPullsPcutstemp->makeTH1D(Form("Z-vertex pull, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"Pull","",BinSettings(200,-10.,10.),Form("h_ZvertPulls_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            //--- The IM and MM histos
            if(i==0 && j==0) h_IMeeg_True = hfOverview->makeTH1D("IM(eeg) True","IM(e^{+}e^{-}#gamma) [MeV]","",BinSettings(500,0.,1000.),"h_IMeeg_True",true);
            h_IMeeg_Rec[i][j] = hfOverviewPcutstemp->makeTH1D(Form("IM(eeg) Rec, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"IM(e^{+}e^{-}#gamma) [MeV]","",BinSettings(500,0.,1000.),Form("h_IMeeg_Rec_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            h_IMeeg_Fit[i][j] = hfOverviewPcutstemp->makeTH1D(Form("IM(eeg) Fit, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"IM(e^{+}e^{-}#gamma) [MeV]","",BinSettings(500,0.,1000.),Form("h_IMeeg_Fit_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            if(i==0 && j==0) h_IMgg_True = hfOverview->makeTH1D("IM(gg) True","IM(#gamma#gamma) [MeV]","",BinSettings(500,0.,1000.),"h_IMgg_True",true);
            h_IMgg_Rec[i][j] = hfOverviewPcutstemp->makeTH1D(Form("IM(gg) Rec, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"IM(#gamma#gamma) [MeV]","",BinSettings(500,0.,1000.),Form("h_IMgg_Rec_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);
            h_IMgg_Fit[i][j] = hfOverviewPcutstemp->makeTH1D(Form("IM(gg) Fit, %s, %s",kftitle[i].c_str(),pcuttitle[j].c_str()),"IM(#gamma#gamma) [MeV]","",BinSettings(500,0.,1000.),Form("h_IMgg_Fit_%s_%s",kfname[i].c_str(),pcutname[j].c_str()),true);

            delete hfFitRecPcutstempCB; delete hfFitRecPcutstempTA; delete hfTrueFitPcutstempCB; delete hfTrueFitPcutstempTA;
            delete hfOverviewPcutstemp; delete hfFitRecPcutstemp; delete hfTrueFitPcutstemp;
        }
        delete hfOverviewKFCasestemp; delete hfFitRecKFCasestemp; delete hfTrueFitKFCasestemp;
    }
}

void scratch_lheijken_checkkinfit::DoMatchTrueRecoStuff(const TParticleList &allmcpart, const std::vector<TParticlePtr> &trueparts, const TCandidateList &recocands, std::vector<std::vector<TParticlePtr>> &matchtruerecpart, double zvert)
{
    //-- One should probably shift the true particle before this matching too
    const auto matched  = utils::match1to1(allmcpart, recocands.get_ptr_list(),
                                           [] (const TParticlePtr& p1, const TCandidatePtr& p2) {return p1->Angle(*p2);},
                                           {0.0, std_ext::degree_to_radian(15.0)});

    vec3 vertshift{0,0,zvert};
    for(auto& truepart: trueparts){

        //-- skip the beam particle
        if(truepart->Type() == ParticleTypeDatabase::BeamTarget) continue;

        TCandidatePtr match = utils::FindMatched(matched,truepart);
        vector<TParticlePtr> truerecpair;
        if(match){
            truerecpair.push_back(truepart);
            int p = -1;
            //-- calculate the true-rec for theta (shifting true with z-vert) and phi
            double kinvar_tr[3]; double kinvar_r[3];
            kinvar_tr[0] = -1000;
            kinvar_tr[1] = ((truepart->p + vertshift).Theta() - match->Theta)*radtodeg;
            kinvar_tr[2] = (truepart->Phi() - match->Phi)*radtodeg;
            //-- and store the reconstructed values of theta and phi
            kinvar_r[0] = -1000; kinvar_r[1] = match->Theta*radtodeg; kinvar_r[2] = match->Phi*radtodeg;
            if(truepart->Type() == ParticleTypeDatabase::Proton){
                truerecpair.push_back(make_shared<TParticle>(ParticleTypeDatabase::Proton,match));
                matchtruerecpart.push_back(truerecpair);
                //-- set the ekin values
                kinvar_tr[0] = truepart->Ek() - truerecpair.at(1)->Ek();
                kinvar_r[0] = truerecpair.at(1)->Ek();
                p = en_p;
            }
            if(truepart->Type() == ParticleTypeDatabase::ePlus){
                truerecpair.push_back(make_shared<TParticle>(ParticleTypeDatabase::ePlus,match));
                matchtruerecpart.push_back(truerecpair);
                //-- set the ekin values
                kinvar_tr[0] = truepart->Ek() - truerecpair.at(1)->Ek();
                kinvar_r[0] = truerecpair.at(1)->Ek();
                p = en_ep;
            }
            if(truepart->Type() == ParticleTypeDatabase::eMinus){
                truerecpair.push_back(make_shared<TParticle>(ParticleTypeDatabase::eMinus,match));
                matchtruerecpart.push_back(truerecpair);
                //-- set the ekin values
                kinvar_tr[0] = truepart->Ek() - truerecpair.at(1)->Ek();
                kinvar_r[0] = truerecpair.at(1)->Ek();
                p = en_em;
            }
            if(truepart->Type() == ParticleTypeDatabase::Photon){
                truerecpair.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,match));
                matchtruerecpart.push_back(truerecpair);
                //-- set the ekin values
                kinvar_tr[0] = truepart->Ek() - truerecpair.at(1)->Ek();
                kinvar_r[0] = truerecpair.at(1)->Ek();
                p = en_g;
            }
            if(p<0) {LOG(INFO)<<"True particle is not p, e+, e- or g"; continue;}
            //-- Fill the true-rec histograms
            if(match->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                for(int i=0; i<nrKinVars;i++){
                    h_TrueRec_vsEk_CB[p][i]->Fill(kinvar_r[0],kinvar_tr[i]);
                    h_TrueRec_vsTh_CB[p][i]->Fill(kinvar_r[1],kinvar_tr[i]);
                    h_TrueRec_vsPh_CB[p][i]->Fill(kinvar_r[2],kinvar_tr[i]);
                }
            }
            else{
                for(int i=0; i<nrKinVars; i++){
                    h_TrueRec_vsEk_TAPS[p][i]->Fill(kinvar_r[0],kinvar_tr[i]);
                    h_TrueRec_vsTh_TAPS[p][i]->Fill(kinvar_r[1],kinvar_tr[i]);
                    h_TrueRec_vsPh_TAPS[p][i]->Fill(kinvar_r[2],kinvar_tr[i]);
                }
            }
        }
    }
}

void scratch_lheijken_checkkinfit::DoFitComparisons(const int kfnr, const int pcnr, const int nrphif, const LorentzVec fitprot, const TParticlePtr recprot, const std::vector<LorentzVec> fitphot, const TParticleList recphot, const std::vector<std::vector<TParticlePtr>> matchtruerecopart, const std::vector<TParticlePtr> truepart, const std::vector<utils::Fitter::FitParticle> fitparts, double zvertfit, double tw)
{
    double kinvar_r[3], kinvar_t[3], kinvar_fr[3], kinvar_tf[3];
    for(int i=0; i<3; i++){kinvar_r[i]=-1000; kinvar_t[i]=-1000; kinvar_fr[i]=-1000; kinvar_tf[i]=-1000;}
    vec3 vertshift{0,0,zvertfit};

    //-- compare proton, if there is a reconstructed proton
    if(recprot){
        //-- calculate fit - rec (for theta: shift the fit with the extracted vertex)
        kinvar_r[0] = recprot->Ek(); kinvar_r[1] = recprot->Theta()*radtodeg; kinvar_r[2] = recprot->Phi()*radtodeg;
        kinvar_fr[0] = (fitprot.E - ParticleTypeDatabase::Proton.Mass()) - kinvar_r[0];
        kinvar_fr[1] = ((fitprot.p + vertshift).Theta()*radtodeg - kinvar_r[1]);
        kinvar_fr[2] = fitprot.Phi()*radtodeg - kinvar_r[2];
        //-- fill the fit-rec histograms
        if(recprot->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            for(int i=0; i<nrKinVars;i++){
                h_FitRec_vsEk_CB[kfnr][pcnr][en_p][i]->Fill(kinvar_r[0],kinvar_fr[i],tw);
                h_FitRec_vsTh_CB[kfnr][pcnr][en_p][i]->Fill(kinvar_r[1],kinvar_fr[i],tw);
                h_FitRec_vsPh_CB[kfnr][pcnr][en_p][i]->Fill(kinvar_r[2],kinvar_fr[i],tw);
            }
        }
        else{
            for(int i=0; i<nrKinVars;i++){
                h_FitRec_vsEk_TAPS[kfnr][pcnr][en_p][i]->Fill(kinvar_r[0],kinvar_fr[i],tw);
                h_FitRec_vsTh_TAPS[kfnr][pcnr][en_p][i]->Fill(kinvar_r[1],kinvar_fr[i],tw);
                h_FitRec_vsPh_TAPS[kfnr][pcnr][en_p][i]->Fill(kinvar_r[2],kinvar_fr[i],tw);
            }
        }
    }
    //-- find the proton among the true particles and calculate true - fit
    TParticlePtr trueprot = 0;
    for(TParticlePtr p: truepart){
        if(p->Type() == ParticleTypeDatabase::Proton){
            trueprot = p;
            break;
        }
    }
    if(trueprot){
        kinvar_t[0] = trueprot->Ek(); kinvar_t[1] = trueprot->Theta()*radtodeg; kinvar_t[2] = trueprot->Phi()*radtodeg;
        kinvar_tf[0] = kinvar_t[0] - (fitprot.E - ParticleTypeDatabase::Proton.Mass());
        kinvar_tf[1] = kinvar_t[1] - fitprot.Theta()*radtodeg;
        kinvar_tf[2] = kinvar_t[2] - fitprot.Phi()*radtodeg;
        //--- fill the true-fit histograms
        if((fitprot.Theta())*radtodeg >= 22.){
            for(int i=0; i<nrKinVars;i++){
                h_TrueFit_vsEk_CB[kfnr][pcnr][en_p][i]->Fill(kinvar_t[0],kinvar_tf[i],tw);
                h_TrueFit_vsTh_CB[kfnr][pcnr][en_p][i]->Fill(kinvar_t[1],kinvar_tf[i],tw);
                h_TrueFit_vsPh_CB[kfnr][pcnr][en_p][i]->Fill(kinvar_t[2],kinvar_tf[i],tw);
            }
        }
        else{
            for(int i=0; i<nrKinVars;i++){
                h_TrueFit_vsEk_TAPS[kfnr][pcnr][en_p][i]->Fill(kinvar_t[0],kinvar_tf[i],tw);
                h_TrueFit_vsTh_TAPS[kfnr][pcnr][en_p][i]->Fill(kinvar_t[1],kinvar_tf[i],tw);
                h_TrueFit_vsPh_TAPS[kfnr][pcnr][en_p][i]->Fill(kinvar_t[2],kinvar_tf[i],tw);
            }
        }
    }
    //-- If the proton was measured, fill the pulls
    int fitpartstart = 0;               // if the proton was unmeasured, the fitparticles only contains the photons
    if(recprot){
        fitpartstart = 1;               // if it was measured, the fitparticles start with the proton and then comes the photons
        if((fitprot.Theta())*radtodeg >= 22.){
            for(int i=0; i<nrFitVars; i++){
                h_PartPulls_CB[kfnr][pcnr][en_p][i]->Fill(fitparts.at(0).GetPulls().at(i),tw);
            }
        }
        else {
            for(int i=0; i<nrFitVars; i++){
                h_PartPulls_TAPS[kfnr][pcnr][en_p][i]->Fill(fitparts.at(0).GetPulls().at(i),tw);
            }
        }
    }

    //-- Do the same with the "photons". If it is MC, classify the rec and fit particles as e+, e-
    //   or g depending on who the rec have been matched with. If it's not, then just fill the g dists
    if(matchtruerecopart.size()>0){
        for(int i=0; i<(int)recphot.size(); i++){
            TParticlePtr rph = recphot.at(i);
            for(int j=0; j<3; j++){kinvar_r[j]=-1000; kinvar_t[j]=-1000; kinvar_fr[j]=-1000; kinvar_tf[j]=-1000;}
            //-- Figure out which particle it might be by comparing the center-element to the true-reco matched
            int ptype=-1;
            TParticlePtr matchrecpart = 0; TParticlePtr matchtruepart = 0;
            for(int j=0; j<(int)matchtruerecopart.size(); j++){
                TParticlePtr tmph = matchtruerecopart.at(j).at(0);
                if(tmph->Type() == ParticleTypeDatabase::Proton) continue;
                TParticlePtr rmph = matchtruerecopart.at(j).at(1);
                if(rph->Candidate->FindCaloCluster()->DetectorType == rmph->Candidate->FindCaloCluster()->DetectorType){
                    if(rph->Candidate->FindCaloCluster()->CentralElement == rmph->Candidate->FindCaloCluster()->CentralElement){
                        matchtruepart = tmph;
                        matchrecpart = rmph;
                    }
                }
            }
            if(matchrecpart){
                if(matchtruepart->Type() == ParticleTypeDatabase::ePlus) ptype = en_ep;
                else if(matchtruepart->Type() == ParticleTypeDatabase::eMinus) ptype = en_em;
                else if(matchtruepart->Type() == ParticleTypeDatabase::Photon) ptype = en_g;
                else {LOG(INFO)<<"The reconstructed photon candidate was not matched up with e+,e- or g";  continue;}

                //-- calculate fit - rec (for theta: shift the fit with the extracted vertex)
                kinvar_r[0] = recphot.at(i)->Ek(); kinvar_r[1] = recphot.at(i)->Theta()*radtodeg; kinvar_r[2] = recphot.at(i)->Phi()*radtodeg;
                kinvar_fr[0] = fitphot.at(i).E - kinvar_r[0];
                kinvar_fr[1] = ((fitphot.at(i).p + vertshift).Theta()*radtodeg - kinvar_r[1]);
                kinvar_fr[2] = fitphot.at(i).Phi()*radtodeg - kinvar_r[2];
                //-- calculate the true - fit
                kinvar_t[0] = matchtruepart->Ek(); kinvar_t[1] = matchtruepart->Theta()*radtodeg; kinvar_t[2] = matchtruepart->Phi()*radtodeg;
                kinvar_tf[0] = kinvar_t[0] - fitphot.at(i).E;
                kinvar_tf[1] = kinvar_t[1] - fitphot.at(i).Theta()*radtodeg;
                kinvar_tf[2] = kinvar_t[2] - fitphot.at(i).Phi()*radtodeg;
                //-- fill the fit-rec histograms
                if(recphot.at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                    for(int j=0; j<nrKinVars;j++){
                        h_FitRec_vsEk_CB[kfnr][pcnr][ptype][j]->Fill(kinvar_r[0],kinvar_fr[j],tw);
                        h_FitRec_vsTh_CB[kfnr][pcnr][ptype][j]->Fill(kinvar_r[1],kinvar_fr[j],tw);
                        h_FitRec_vsPh_CB[kfnr][pcnr][ptype][j]->Fill(kinvar_r[2],kinvar_fr[j],tw);
                    }
                }
                else{
                    for(int j=0; j<nrKinVars;j++){
                        h_FitRec_vsEk_TAPS[kfnr][pcnr][ptype][j]->Fill(kinvar_r[0],kinvar_fr[j],tw);
                        h_FitRec_vsTh_TAPS[kfnr][pcnr][ptype][j]->Fill(kinvar_r[1],kinvar_fr[j],tw);
                        h_FitRec_vsPh_TAPS[kfnr][pcnr][ptype][j]->Fill(kinvar_r[2],kinvar_fr[j],tw);
                    }
                }
                //-- fill true-fit histograms
                if((fitphot.at(i).Theta())*radtodeg >= 22.0){
                    for(int j=0; j<nrKinVars;j++){
                        h_TrueFit_vsEk_CB[kfnr][pcnr][ptype][j]->Fill(kinvar_t[0],kinvar_tf[j],tw);
                        h_TrueFit_vsTh_CB[kfnr][pcnr][ptype][j]->Fill(kinvar_t[1],kinvar_tf[j],tw);
                        h_TrueFit_vsPh_CB[kfnr][pcnr][ptype][j]->Fill(kinvar_t[2],kinvar_tf[j],tw);
                    }
                }
                else{
                    for(int j=0; j<nrKinVars;j++){
                        h_TrueFit_vsEk_TAPS[kfnr][pcnr][ptype][j]->Fill(kinvar_t[0],kinvar_tf[j],tw);
                        h_TrueFit_vsTh_TAPS[kfnr][pcnr][ptype][j]->Fill(kinvar_t[1],kinvar_tf[j],tw);
                        h_TrueFit_vsPh_TAPS[kfnr][pcnr][ptype][j]->Fill(kinvar_t[2],kinvar_tf[j],tw);
                    }
                }
                //-- fill the pull histograms
                if((fitphot.at(i).Theta())*radtodeg >= 22.){
                    for(int j=0; j<nrFitVars; j++){
                        h_PartPulls_CB[kfnr][pcnr][ptype][j]->Fill(fitparts.at(fitpartstart+i).GetPulls().at(j),tw);
                    }
                }
                else {
                    for(int j=0; j<nrFitVars; j++){
                        h_PartPulls_TAPS[kfnr][pcnr][ptype][j]->Fill(fitparts.at(fitpartstart+i).GetPulls().at(j),tw);
                    }
                }
            }
        }
    }
    //-- If it's not MC then just fill the photon distribution (could try to think of something clever)
    //-- calculate fit - rec (for theta: shift the fit with the extracted vertex)
    else{
        for(int i=0; i<(int)recphot.size(); i++){
            kinvar_r[0] = recphot.at(i)->Ek(); kinvar_r[1] = recphot.at(i)->Theta()*radtodeg; kinvar_r[2] = recphot.at(i)->Phi()*radtodeg;
            kinvar_fr[0] = fitphot.at(i).E - kinvar_r[0];
            kinvar_fr[1] = ((fitphot.at(i).p + vertshift).Theta()*radtodeg - kinvar_r[1]);
            kinvar_fr[2] = fitphot.at(i).Phi()*radtodeg - kinvar_r[2];
            //-- fill the fit-rec histograms
            if(recphot.at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                for(int j=0; j<nrKinVars;j++){
                    h_FitRec_vsEk_CB[kfnr][pcnr][en_g][j]->Fill(kinvar_r[0],kinvar_fr[j],tw);
                    h_FitRec_vsTh_CB[kfnr][pcnr][en_g][j]->Fill(kinvar_r[1],kinvar_fr[j],tw);
                    h_FitRec_vsPh_CB[kfnr][pcnr][en_g][j]->Fill(kinvar_r[2],kinvar_fr[j],tw);
                }
            }
            else{
                for(int j=0; j<nrKinVars;j++){
                    h_FitRec_vsEk_TAPS[kfnr][pcnr][en_g][j]->Fill(kinvar_r[0],kinvar_fr[j],tw);
                    h_FitRec_vsTh_TAPS[kfnr][pcnr][en_g][j]->Fill(kinvar_r[1],kinvar_fr[j],tw);
                    h_FitRec_vsPh_TAPS[kfnr][pcnr][en_g][j]->Fill(kinvar_r[2],kinvar_fr[j],tw);
                }
            }
            //-- fill the pull histograms
            if(recphot.at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                for(int j=0; j<nrFitVars; j++){
                    h_PartPulls_CB[kfnr][pcnr][en_g][j]->Fill(fitparts.at(fitpartstart+i).GetPulls().at(j),tw);
                }
            }
            else {
                for(int j=0; j<nrFitVars; j++){
                    h_PartPulls_TAPS[kfnr][pcnr][en_g][j]->Fill(fitparts.at(fitpartstart+i).GetPulls().at(j),tw);
                }
            }
        }
    }

    //-- Fill the IM histograms
    if(nrphif == 3){
        h_IMeeg_Fit[kfnr][pcnr]->Fill((fitphot.at(0) + fitphot.at(1) + fitphot.at(2)).M(),tw);
        h_IMeeg_Rec[kfnr][pcnr]->Fill((*recphot.at(0) + *recphot.at(1) + *recphot.at(2)).M(),tw);
    }
    if(nrphif == 2){
        h_IMgg_Fit[kfnr][pcnr]->Fill((fitphot.at(0) + fitphot.at(1)).M(),tw);
        h_IMgg_Rec[kfnr][pcnr]->Fill((*recphot.at(0) + *recphot.at(1)).M(),tw);
    }
    if(nrphif == 4){
        for(int i=0; i<3; i++){
            for(int j=(i+1);j<4; j++){
                h_IMgg_Fit[kfnr][pcnr]->Fill((fitphot.at(i) + fitphot.at(j)).M(),tw);
                h_IMgg_Rec[kfnr][pcnr]->Fill((*recphot.at(i) + *recphot.at(j)).M(),tw);
            }
        }
    }
}

APLCON::Fit_Settings_t scratch_lheijken_checkkinfit::MakeFitSettings(unsigned max_iterations)
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    //settings.ConstraintAccuracy = 1.0e-10;
    //settings.Chi2Accuracy = 1.0e-8;
    //settings.DebugLevel = 0;
    return settings;
}

AUTO_REGISTER_PHYSICS(scratch_lheijken_checkkinfit)
