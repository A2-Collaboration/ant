#include "Pi0Dalitz.h"

#include "base/ParticleType.h"
#include "plot/HistogramFactory.h"
#include "utils/Combinatorics.h"
#include "utils/ParticleTools.h"
#include "utils/Matcher.h"
#include "base/ParticleTypeTree.h"
#include "utils/uncertainties/Interpolated.h"
#include "utils/uncertainties/FitterSergey.h"
#include "expconfig/ExpConfig.h"
#include "base/Logger.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "base/WrapTFile.h"

#include <iostream>
#include <memory>

/**
 * Analysis of the pi0->e+e-g channel
 *
 **/

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

Pi0Dalitz::Pi0Dalitz(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    fit_model(utils::UncertaintyModels::Interpolated::makeAndLoad(
                             utils::UncertaintyModels::Interpolated::Type_t::MC,
                             make_shared<utils::UncertaintyModels::FitterSergey>())),
    fitter(nullptr, opts->Get<bool>("FitZVertex", true)),
    npfitter(nullptr, opts->Get<bool>("FitZVertex", true),MakeFitSettings(30))
{
    fitter.SetZVertexSigma(3.0);
    npfitter.SetZVertexSigma(3.0);

    tagger_detector = ExpConfig::Setup::GetDetector<expconfig::detector::Tagger>();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    pid_detector = ExpConfig::Setup::GetDetector<expconfig::detector::PID>();
    veto_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPSVeto>();

    CreateHistos();
    TFFTree.CreateBranches(HistFac.makeTTree("tff_variables"));

    promptrandom.AddPromptRange({-8, 8});
    promptrandom.AddRandomRange({-25, -15});
    promptrandom.AddRandomRange({15, 25});

    //-- Create the names of the files containing the cluster-E corrections and load the respective correction histograms
    if(doecorr){
        string lepCBfile = std_ext::formatter() << ExpConfig::Setup::Get().GetPhysicsFilesDirectory() << "/Corr_h_Ek_TrueRecvsRec_em_CB.root";
        lepCBcorr.LoadECorr(lepCBfile,"hCorrectionsFinal");
        string photCBfile = std_ext::formatter() << ExpConfig::Setup::Get().GetPhysicsFilesDirectory() << "/Corr_h_Ek_TrueRecvsRec_g_CB.root";
        photCBcorr.LoadECorr(photCBfile,"hCorrectionsFinal");
        string photTAPSfile = std_ext::formatter() << ExpConfig::Setup::Get().GetPhysicsFilesDirectory() << "/Corr_h_Ek_TrueRecvsRec_g_TAPS.root";
        photTAPScorr.LoadECorr(photTAPSfile,"hCorrectionsFinal");
    }

}


void Pi0Dalitz::ProcessEvent(const TEvent& event, manager_t&)
{
    triggersimu.ProcessEvent(event);

    //-- Check the decay string for MC pattern
    bool isMC=false;
    string decay;
    if(event.Reconstructed().ID.isSet(ant::TID::Flags_t::MC) && event.MCTrue().ParticleTree){
        decay = utils::ParticleTools::GetDecayString(event.MCTrue().ParticleTree);
        isMC=true;
    }
    else
        decay = "data" + ExpConfig::Setup::Get().GetName();
    vector<bool> WhichMC; // 0 : pi0->eeg, 1 : pi0->gg, 2 : all other MC
    bool MCpi0eeg = false; bool MCpi0gg = false; bool MCother = false;
    if(isMC){
        MCpi0eeg  = event.MCTrue().ParticleTree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_eeg),
                                                         utils::ParticleTools::MatchByParticleName);
        MCpi0gg  = event.MCTrue().ParticleTree->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g),
                                                        utils::ParticleTools::MatchByParticleName);
        if(!(MCpi0eeg || MCpi0gg))
            MCother = true;
    }
    WhichMC.push_back(MCpi0eeg); WhichMC.push_back(MCpi0gg); WhichMC.push_back(MCother);

    //-- The reconstructed candidates
    const auto& candidates = event.Reconstructed().Candidates;
    h_AnalysisStat->Fill(1.);
    //-- Make a selection of the candidates and sort them into charged and neutral and apply energy corrections
    //   Also check if they fullfil the requirements of good lepton/photon/proton candidates and if the event
    //   fullfil the requirement to be a Dalitz decay or an 2 photon decay event
    TCandidatePtrList selectedcands;
    vector<TParticlePtr> protons_to_fit;
    vector<TParticleList> photons_to_fit;
    vector<bool> cases;
    DoCandSelStuff(candidates,selectedcands,protons_to_fit,photons_to_fit,cases);
    bool HasProtonCand = cases.at(0);
    bool IsDalDecay = cases.at(1);
    bool Is2gDecay = cases.at(2);

    //-- Set fit models, np = no proton
    fitter.SetUncertaintyModel(fit_model);
    npfitter.SetUncertaintyModel(fit_model);

    //-- MC stuff, for the pi0 -> eeg and pi0->gg channels
    vector<TParticlePtr> TruePart;
    vector<TParticlePtr> RecoMatchPart;
    if(isMC){
        for(int i=0; i<nrPartTypes; i++) h_AnalysisStat_RecMat->Fill(0.,(double)i);
        //--- Fetch the true final state particles
        if(MCpi0eeg){
            auto TruePartp_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton,event.MCTrue().ParticleTree);
            auto TruePartep_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::ePlus,event.MCTrue().ParticleTree);
            auto TruePartem_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::eMinus,event.MCTrue().ParticleTree);
            auto TruePartg_eeg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Photon,event.MCTrue().ParticleTree);
            TruePart.push_back(TruePartp_eeg); TruePart.push_back(TruePartep_eeg);
            TruePart.push_back(TruePartem_eeg); TruePart.push_back(TruePartg_eeg);
        }
        if(MCpi0gg){
            auto TruePartp_gg = utils::ParticleTools::FindParticle(ParticleTypeDatabase::Proton,event.MCTrue().ParticleTree);
            auto TruePartgs_gg = utils::ParticleTools::FindParticles(ParticleTypeDatabase::Photon, event.MCTrue().ParticleTree);
            if(TruePartgs_gg.size() != 2){LOG(INFO)<<"Didn't find two photons in pi0->gg trueMC"; return;}
            TruePart.push_back(TruePartp_gg);
            TruePart.push_back(TruePartgs_gg.front()); TruePart.push_back(TruePartgs_gg.back());
        }

        //--- Match the true particles to the reconstructed candidates
        auto mcparticleslist = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);
        const auto mcparticles = mcparticleslist.GetAll();
        DoMatchTrueRecoStuff(mcparticles, TruePart, selectedcands, RecoMatchPart);

        //--- For any other MC channel, simply store all produced particles in the truepart vector
        auto alltrueparts = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree).GetAll();
        for(auto& truepart : alltrueparts){
            if(truepart->Type() == ParticleTypeDatabase::BeamTarget) continue;
            TruePart.push_back(truepart);
        }
    }

    //-- Loop over the tagger hits

    //--- The result from the resulting kinfit of each tagger hit, to be stored in the tree
    vector<double> kfimeegs;
    vector<double> kfimees;
    vector<double> kfimggs;
    vector<double> kfprobs;
    vector<bool> eesamepids;
    vector<double> taggerweights;

    for (const auto &tc : event.Reconstructed().TaggerHits) {
        double cortagtime = triggersimu.GetCorrectedTaggerTime(tc);
        promptrandom.SetTaggerTime(cortagtime);
        if (promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        double taggweight = promptrandom.FillWeight();
        TLorentzVector InitialPhotonVec = tc.GetPhotonBeam();

        //-- No cuts (except the intrinsic requirement of a prompt taggerhit)
        DoTaggerStuff(en_nocut,InitialPhotonVec,tc.Time,cortagtime,taggweight);
        DoTriggerStuff(en_nocut,taggweight);
        DoRecoCandStuff(en_nocut,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
        DoTrueMCStuff(en_nocut,WhichMC,TruePart,taggweight);

        //-- CBEsum cut (only on MC, data is always true)
        if(!triggersimu.HasTriggered()) continue;
        DoTaggerStuff(en_cbsum,InitialPhotonVec,tc.Time,cortagtime,taggweight);
        DoTriggerStuff(en_cbsum,taggweight);
        DoRecoCandStuff(en_cbsum,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
        DoTrueMCStuff(en_cbsum,WhichMC,TruePart,taggweight);

        //-- Only allow tagger hits with E < 450 MeV
        if(InitialPhotonVec.E() > 450.) continue;
        DoTaggerStuff(en_tagge,InitialPhotonVec,tc.Time,cortagtime,taggweight);
        DoTriggerStuff(en_tagge,taggweight);
        DoRecoCandStuff(en_tagge,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
        DoTrueMCStuff(en_tagge,WhichMC,TruePart,taggweight);

        //-- Fit the Dalitz decay case with a proton candidate
        if(IsDalDecay && HasProtonCand){
            DoTaggerStuff(en_1p3g,InitialPhotonVec,tc.Time,cortagtime,taggweight);
            DoTriggerStuff(en_1p3g,taggweight);
            DoRecoCandStuff(en_1p3g,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
            DoTrueMCStuff(en_1p3g,WhichMC,TruePart,taggweight);
            //-- Do the kinfit on the hypothesis of measured proton (except its energy) and 3 "photons"
            vector<TParticlePtr> bestprobrec;
            vector<TParticlePtr> bestprobfit;
            double kfprob = DoKinFitStuffWProt(0,protons_to_fit,photons_to_fit,InitialPhotonVec,fitter,bestprobrec,bestprobfit,taggweight);
            //-- Fill all check-histograms if the kf-probability is higher than a given value
            if(kfprob > evselhistoprobcut){
                DoTaggerStuff(en_1p3g_pc,InitialPhotonVec,tc.Time,cortagtime,taggweight);
                DoTriggerStuff(en_1p3g_pc,taggweight);
                DoRecoCandStuff(en_1p3g_pc,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
                DoTrueMCStuff(en_1p3g_pc,WhichMC,TruePart,taggweight);
            }
            //-- If it passes the kinfit, do a fit without protons on the selected photons
            if(kfprob > -50){
                //-- Prepare the photons selected in the previous fit to a format expected by the next fit
                vector<TParticleList> photons_to_fit_2;
                TParticleList phtemp;
                for(const auto &ph : bestprobrec){
                    if(ph->Type() == ParticleTypeDatabase::Photon)
                        phtemp.push_back(ph);
                }
                photons_to_fit_2.push_back(phtemp);
                if(photons_to_fit_2.at(0).size() == 3){
                    //-- Do the kinfit on the hypothesis of an unmeasured proton and 3 "photons"
                    TParticleList fitted_photons;
                    LorentzVec fitted_proton;
                    double kfprobnp = DoKinFitStuffNProt(1,photons_to_fit_2,InitialPhotonVec,npfitter,fitted_photons,fitted_proton,taggweight);
                    //-- Fill all check-histograms if the kf-probability is higher than a given value
                    if(kfprobnp > evselhistoprobcut){
                        DoTaggerStuff(en_1p3gnp_pc,InitialPhotonVec,tc.Time,cortagtime,taggweight);
                        DoTriggerStuff(en_1p3gnp_pc,taggweight);
                        DoRecoCandStuff(en_1p3gnp_pc,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
                        DoTrueMCStuff(en_1p3gnp_pc,WhichMC,TruePart,taggweight);
                    }
                    //-- If it passes the kin fit, select the e+e- pair and store variables for the tree
                    if(kfprobnp > -50){
                        TParticlePtr kfphoton, kfe1, kfe2;
                        bool firste = true;
                        for(const auto& kfpart : fitted_photons){
                            if(kfpart->Type() == ParticleTypeDatabase::Photon){
                                if(kfpart->Candidate->VetoEnergy <= pidecut_chaneu)
                                    kfphoton = kfpart;
                                else {
                                    if(firste) {
                                        kfe1 = kfpart;
                                        firste = false;
                                    }
                                    else
                                        kfe2 = kfpart;
                                }
                            }
                        }
                        //-- Check if the chosen ee pair has the same PID element
                        bool eesamepid = false;
                        if(kfe1->Candidate->FindVetoCluster()->CentralElement == kfe2->Candidate->FindVetoCluster()->CentralElement)
                            eesamepid = true;

                        kfimees.push_back((*kfe1 + *kfe2).M());
                        kfimeegs.push_back((*kfphoton + *kfe1 + *kfe2).M());
                        eesamepids.push_back(eesamepid);
                        kfprobs.push_back(kfprobnp);
                        taggerweights.push_back(taggweight);
                    }
                }
                else
                    LOG(INFO)<<"There should be 3 photons to fit for the Dalitz decay case";
            }
        }

        //-- Fit the Dalitz decay case with no reconstructed proton candidate
        if(IsDalDecay && !HasProtonCand){
            DoTaggerStuff(en_np3g,InitialPhotonVec,tc.Time,cortagtime,taggweight);
            DoTriggerStuff(en_np3g,taggweight);
            DoRecoCandStuff(en_np3g,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
            DoTrueMCStuff(en_np3g,WhichMC,TruePart,taggweight);
            //-- Do the kinfit on the hypothesis of an unmeasured proton and 3 "photons"
            TParticleList fitted_photons;
            LorentzVec fitted_proton;
            double kfprobnp = DoKinFitStuffNProt(2,photons_to_fit,InitialPhotonVec,npfitter,fitted_photons,fitted_proton,taggweight);
            //-- Fill all check-histograms if the kf-probability is higher than a given value
            if(kfprobnp > evselhistoprobcut){
                DoTaggerStuff(en_np3g_pc,InitialPhotonVec,tc.Time,cortagtime,taggweight);
                DoTriggerStuff(en_np3g_pc,taggweight);
                DoRecoCandStuff(en_np3g_pc,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
                DoTrueMCStuff(en_np3g_pc,WhichMC,TruePart,taggweight);
            }
            //-- If it passes the kin fit, select the e+e- pair and store variables for the tree
            if(kfprobnp > -50){
                TParticlePtr kfphoton, kfe1, kfe2;
                bool firste = true;
                for(const auto& kfpart : fitted_photons){
                    if(kfpart->Candidate->VetoEnergy <= pidecut_chaneu)
                        kfphoton = kfpart;
                    else {
                        if(firste) {
                            kfe1 = kfpart;
                            firste = false;
                        }
                        else
                            kfe2 = kfpart;
                    }
                }
                //-- Check if the chosen ee pair has the same PID element
                bool eesamepid = false;
                if(kfe1->Candidate->FindVetoCluster()->CentralElement == kfe2->Candidate->FindVetoCluster()->CentralElement)
                    eesamepid = true;

                kfimees.push_back((*kfe1 + *kfe2).M());
                kfimeegs.push_back((*kfphoton + *kfe1 + *kfe2).M());
                eesamepids.push_back(eesamepid);
                kfprobs.push_back(kfprobnp);
                taggerweights.push_back(taggweight);
            }
        }

        //-- Fit the 2gamma decay case with a proton candidate
        if(Is2gDecay && HasProtonCand){
            DoTaggerStuff(en_1p2g,InitialPhotonVec,tc.Time,cortagtime,taggweight);
            DoTriggerStuff(en_1p2g,taggweight);
            DoRecoCandStuff(en_1p2g,selectedcands,protons_to_fit,photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
            DoTrueMCStuff(en_1p2g,WhichMC,TruePart,taggweight);
            //-- Do the kinfit on the hypothesis of measured proton (except its energy) and 2 photons
            vector<TParticlePtr> bestprobrec;
            vector<TParticlePtr> bestprobfit;
            double kfprob = DoKinFitStuffWProt(3,protons_to_fit,photons_to_fit,InitialPhotonVec,fitter,bestprobrec,bestprobfit,taggweight);
            //-- Fill all check-histograms if the kf-probability is higher than a given value
            if(kfprob > evselhistoprobcut){
                DoTaggerStuff(en_1p2g_pc,InitialPhotonVec,tc.Time,cortagtime,taggweight);
                DoTriggerStuff(en_1p2g_pc,taggweight);
                DoRecoCandStuff(en_1p2g_pc,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
                DoTrueMCStuff(en_1p2g_pc,WhichMC,TruePart,taggweight);
            }
            //-- If it passes the kinfit, do a fit without protons on the selected photons
            if(kfprob > -50){
                //-- Prepare the photons selected in the previous fit to a format expected by the next fit
                vector<TParticleList> photons_to_fit_2;
                TParticleList phtemp;
                for(const auto &ph : bestprobrec){
                    if(ph->Type() == ParticleTypeDatabase::Photon)
                        phtemp.push_back(ph);
                }
                photons_to_fit_2.push_back(phtemp);
                if(photons_to_fit_2.at(0).size() == 2){
                    //-- Do the kinfit on the hypothesis of an unmeasured proton and 2 "photons"
                    TParticleList fitted_photons;
                    LorentzVec fitted_proton;
                    double kfprobnp = DoKinFitStuffNProt(4,photons_to_fit_2,InitialPhotonVec,npfitter,fitted_photons,fitted_proton,taggweight);
                    //-- Fill all check-histograms if the kf-probability is higher than a given value
                    if(kfprobnp > evselhistoprobcut){
                        DoTaggerStuff(en_1p2gnp_pc,InitialPhotonVec,tc.Time,cortagtime,taggweight);
                        DoTriggerStuff(en_1p2gnp_pc,taggweight);
                        DoRecoCandStuff(en_1p2gnp_pc,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
                        DoTrueMCStuff(en_1p2gnp_pc,WhichMC,TruePart,taggweight);
                    }
                    //-- If it passes the kin fit, store variables for the tree
                    if(kfprobnp > -50){
                        kfimggs.push_back((*fitted_photons.at(0) + *fitted_photons.at(1)).M());
                        kfprobs.push_back(kfprobnp);
                        taggerweights.push_back(taggweight);
                    }
                }
                else
                    LOG(INFO)<<"There should be 2 photons to fit for the 2gamma decay case";
            }
        }

        //-- Fit the 2g decay case with no reconstructed proton candidate
        if(Is2gDecay && !HasProtonCand){
            DoTaggerStuff(en_np2g,InitialPhotonVec,tc.Time,cortagtime,taggweight);
            DoTriggerStuff(en_np2g,taggweight);
            DoRecoCandStuff(en_np2g,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
            DoTrueMCStuff(en_np2g,WhichMC,TruePart,taggweight);
            //-- Do the kinfit on the hypothesis of an unmeasured proton and 3 "photons"
            TParticleList fitted_photons;
            LorentzVec fitted_proton;
            double kfprobnp = DoKinFitStuffNProt(5,photons_to_fit,InitialPhotonVec,npfitter,fitted_photons,fitted_proton,taggweight);
            //-- Fill all check-histograms if the kf-probability is higher than a given value
            if(kfprobnp > evselhistoprobcut){
                DoTaggerStuff(en_np2g_pc,InitialPhotonVec,tc.Time,cortagtime,taggweight);
                DoTriggerStuff(en_np2g_pc,taggweight);
                DoRecoCandStuff(en_np2g_pc,selectedcands,protons_to_fit, photons_to_fit,RecoMatchPart,WhichMC,InitialPhotonVec,taggweight);
                DoTrueMCStuff(en_np2g_pc,WhichMC,TruePart,taggweight);
            }
            //-- If it passes the kin fit, store variables for the tree
            if(kfprobnp > -50){
                kfimggs.push_back((*fitted_photons.at(0) + *fitted_photons.at(1)).M());
                kfprobs.push_back(kfprobnp);
                taggerweights.push_back(taggweight);
            }
        }
    }

    //--- Fill the TFF tree
    if((int)kfimees.size()>0 || (int)kfimggs.size()>0){
        TFFTree.TBIsDalDec = IsDalDecay;
        TFFTree.TBIs2gDec = Is2gDecay;
        TFFTree.TBHasProtonCand = HasProtonCand;
        TFFTree.TBIMees = kfimees;
        TFFTree.TBIMeegs = kfimeegs;
        TFFTree.TBIMggs = kfimggs;
        TFFTree.TBsamePIDs = eesamepids;
        TFFTree.TBKFprobs = kfprobs;
        TFFTree.TBTaggWeights = taggerweights;
        TFFTree.Tree->Fill();
    }
}

void Pi0Dalitz::Finish()
{
    //-- After-treatment of statistics histograms
    double stat = h_RecoTrueMatch->GetBinContent(1);
    double eff;
    if(stat > 0){
        for(int i=1; i<= h_RecoTrueMatch->GetNbinsX(); i++){
            eff = h_RecoTrueMatch->GetBinContent(i)/stat;
            h_RecoTrueMatch->SetBinContent(i,eff);
            h_RecoTrueMatch->SetBinError(i,0); //until I find a better way of dealing with the errors
        }
    }
}

void Pi0Dalitz::ShowResult()
{

}

void Pi0Dalitz::CreateHistos()
{
    //-- Histogram folders
    auto hfTrueMC = new HistogramFactory("TrueMC",HistFac,"");
    auto hfTrigChecks = new HistogramFactory("TrigChecks",HistFac,"");
    auto hfTaggChecks = new HistogramFactory("TaggChecks",HistFac,"");
    auto hfCandChecks = new HistogramFactory("CandChecks",HistFac,"");
    auto hfOverview = new HistogramFactory("Overview",HistFac,"");
    auto hfKinFit = new HistogramFactory("KinFit",HistFac,"");

    //-- Specific bin settings
    //--- Tagger: make some bin edges half way between the energy values.
    vector<double> photonenergies;
    int ntaggch = tagger_detector->GetNChannels();
    for(int ch=0;ch<ntaggch;ch++)
        photonenergies.push_back(tagger_detector->GetPhotonEnergy(ch));
    vector<double> energy_bins(ntaggch+1);
    for(int b=0;b<(ntaggch-1);b++){
      energy_bins.at(b+1)=0.5*(photonenergies.at(b)+photonenergies.at(b+1));
    }
    //--- Tagger: and the top and bottom have width of the adjacent bin
    energy_bins.at(0) = energy_bins.at(1) - (energy_bins.at(2) - energy_bins.at(1));
    energy_bins.at(ntaggch) = energy_bins.at(ntaggch-1) + (energy_bins.at(ntaggch-1)- energy_bins.at(ntaggch-2));

    int npidch = pid_detector->GetNChannels();
    int nvetoch = veto_detector->GetNChannels();

    //-- The histograms
    //--- MC true
    h_IMeegTrue = hfTrueMC->makeTH1D("IMeeg True","IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),"h_IMeegTrue",true);
    h_IMggTrue = hfTrueMC->makeTH1D("IMgg True","IM(#gamma#gamma)","",BinSettings(250,0.,1000.),"h_IMggTrue",true);
    h_IMeeTrue = hfTrueMC->makeTH1D("IMee True","IM(e^{+}e^{-})","",BinSettings(150,0.,150.),"h_IMeeTrue",true);
    h_RecoTrueMatch = hfTrueMC->makeTH1D("Reco-True matches","Particle","Fraction where a matched reco candidate was found",BinSettings(5,-0.5,4.5),"h_RecoTrueMatch",true);
    h_RecoTrueAngle = hfTrueMC->makeTH2D("Angle between matched reco cand and true mc","Particle","opening angle [deg]",BinSettings(5,-0.5,4.5),BinSettings(180,-0.5,44.5),"h_RecoTrueAngle",true);
    string particlename[] = {"p","eplus","eminus","photon"};
    string particletitle[] = {"p","e^{+}","e^{-}","#gamma"};
    h_RecoTrueMatch->GetXaxis()->SetBinLabel(1,"# events");
    for(int i=0; i<nrPartTypes; i++){
        h_RecoTrueMatch->GetXaxis()->SetBinLabel(i+2,particlename[i].c_str());
        h_RecoTrueAngle->GetXaxis()->SetBinLabel(i+2,particlename[i].c_str());
        h_EktrueEkrec[i] = hfTrueMC->makeTH1D(Form("Ek_{true}/Ek_{rec} %s",particlename[i].c_str()),"Ek_{true}/Ek_{rec}","",BinSettings(100,0.5,1.5),Form("h_EktrueEkrec%s",particlename[i].c_str()),true);
    }
    h_EktrueEkrec_gg = hfTrueMC->makeTH1D("Ek_{true}/Ek_{rec} for #gamma pairs","Ek_{true}/Ek_{rec}","",BinSettings(100,0.5,1.5),"h_EktrueEkrec_gg",true);

    //--- Reconstructed candidates info from candidate selection
    auto hfCandSelect = new HistogramFactory("CandSel",*hfCandChecks,"");
    string candselname[] = {"All","Time","CBNeu","CBCha","TAPSNeu","TAPSCha","DalDec","DalDecPhotLep","DalDecProt","TwoPhotDec","TwoPhotDecPhot","TwoPhotDecProt"};
    string candseltitle[] = {"All","Time","CB - Neutral","CB - Charged","TAPS - Neutral","TAPS - Charged","Dalitz decay","Dalitz decay - photons/leptons","Dalitz decay - protons","Two photon decay","Two photon decay - photons","Two photon decay - protons"};
    vector<HistogramFactory> hfCandSels;
    for(int i=0; i<nrCandSel; i++){
        //---- Create subfolder for each candidate selection
        auto hfCandSelsTemp = new HistogramFactory(candselname[i],*hfCandSelect,"");
        hfCandSels.push_back(*hfCandSelsTemp);
        //---- Histos
        h_PIDEvsT_cs[i] = hfCandSels.at(i).makeTH2D(Form("PID E vs T, %s",candseltitle[i].c_str()),"PID Energy","PID Time",BinSettings(200,0.,20),BinSettings(241,-60.5,60.5),Form("h_PIDEvsT%s",candselname[i].c_str()),true);
        h_TVetoEvsT_cs[i] = hfCandSels.at(i).makeTH2D(Form("TAPSVeto E vs T, %s",candseltitle[i].c_str()),"TAPSVeto Energy","TAPSVeto Time",BinSettings(200,0.,20),BinSettings(241,-60.5,60.5),Form("h_TVetoEvsT%s",candselname[i].c_str()),true);
        h_CBEvsT_cs[i] = hfCandSels.at(i).makeTH2D(Form("CB E vs T, %s",candseltitle[i].c_str()),"CB Energy","CB Time",BinSettings(200,0.,600),BinSettings(241,-60.5,60.5),Form("h_CBEvsT%s",candselname[i].c_str()),true);
        h_TAPSEvsT_cs[i] = hfCandSels.at(i).makeTH2D(Form("TAPS E vs T, %s",candseltitle[i].c_str()),"TAPS Energy","TAPS Time",BinSettings(200,0.,600),BinSettings(241,-60.5,60.5),Form("h_TAPSEvsT%s",candselname[i].c_str()),true);
        h_CBEvsNrCr_cs[i] = hfCandSels.at(i).makeTH2D(Form("CB E vs NrCr, %s",candseltitle[i].c_str()),"CB Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_CBEvsNrCr%s",candselname[i].c_str()),true);
        h_TAPSEvsNrCr_cs[i] = hfCandSels.at(i).makeTH2D(Form("TAPS E vs NrCr, %s",candseltitle[i].c_str()),"TAPS Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_TAPSEvsNrCr%s",candselname[i].c_str()),true);
        h_TimeCBvsPID_cs[i] = hfCandSels.at(i).makeTH2D(Form("Time CB vs PID, %s",candseltitle[i].c_str()),"CB Time","PID Time",BinSettings(241,-60.5,60.5),BinSettings(241,-60.5,60.5),Form("h_TimeCBvsPID%s",candselname[i].c_str()),true);
        h_EnergyCBvsPID_cs[i] = hfCandSels.at(i).makeTH2D(Form("Energy CB vs PID, %s",candseltitle[i].c_str()),"CB Energy","PID Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyCBvsPID%s",candselname[i].c_str()),true);
        h_TimeTAPSvsTVeto_cs[i] = hfCandSels.at(i).makeTH2D(Form("Time TAPS vs TAPSVeto, %s",candseltitle[i].c_str()),"TAPS Time","TAPSVeto Time",BinSettings(241,-60.5,60.5),BinSettings(241,-60.5,60.5),Form("h_TimeTAPSvsTAPSVeto%s",candselname[i].c_str()),true);
        h_EnergyTAPSvsTVeto_cs[i] = hfCandSels.at(i).makeTH2D(Form("Energy TAPS vs TAPSVeto, %s",candseltitle[i].c_str()),"TAPS Energy","TAPSVeto Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyTAPSvsTAPSVeto%s",candselname[i].c_str()),true);
    }

    //--- Reconstructed candidates info from event selection
    h_PIDMultUsed = hfCandChecks->makeTH2D("Multiple usage of a PID element","Nr times used/event","Element nr",BinSettings(10,-0.5,9.5),BinSettings(npidch),"h_PIDMultUsed",true);
    h_VetoMultUsed = hfCandChecks->makeTH2D("Multiple usage of a Veto element","Nr times used/event","Element nr",BinSettings(10,-0.5,9.5),BinSettings(nvetoch),"h_VetoMultUsed",true);

    string cutname[] = {"NoCuts","CBE","Eg","DalDecwProt","DalDecwProtPcut","DalDecwProtNPfitPcut","DalDecnProt","DalDecnProtPcut","2gDecwProt","2gDecwProtPcut","2gDecwProtNPfitPcut","2gDecnProt","2gDecnProtPcut"};
    string cuttitle[] = {"no cuts","CBEsum","Eg","Dalitz decay with proton","Dalitz decay with proton, P(#chi^{2})>"+std::to_string(evselhistoprobcut),"Dalitz decay with proton, noprot-fit, P(#chi^{2})>"+std::to_string(evselhistoprobcut),"Dalitz decay no proton","Dalitz decay no proton, P(#chi^{2})>"+std::to_string(evselhistoprobcut),"2#gamma with proton","2#gamma with proton, P(#chi^{2})>"+std::to_string(evselhistoprobcut),"2#gamma with proton, noprot-fit, P(#chi^{2})>"+std::to_string(evselhistoprobcut),"2#gamma no proton","2#gamma no proton, P(#chi^{2})>"+std::to_string(evselhistoprobcut)};
    vector<HistogramFactory> hfTrigChecksCuts, hfTaggChecksCuts, hfCandChecksCuts, hfTrueChecksCuts;
    for(int i=0; i<nrEvSel; i++){
        //---- Create subfolders for each cut
        auto hfTrigChecksCutTemp = new HistogramFactory(cutname[i],*hfTrigChecks,"");
        hfTrigChecksCuts.push_back(*hfTrigChecksCutTemp);
        auto hfTaggChecksCutTemp = new HistogramFactory(cutname[i],*hfTaggChecks,"");
        hfTaggChecksCuts.push_back(*hfTaggChecksCutTemp);
        auto hfCandChecksCutTemp = new HistogramFactory(cutname[i],*hfCandChecks,"");
        hfCandChecksCuts.push_back(*hfCandChecksCutTemp);
        auto hfMatchedTemp = new HistogramFactory("ForMCMatched",hfCandChecksCuts.at(i),"");
        auto hfTrueChecksCutTemp = new HistogramFactory(cutname[i],*hfTrueMC,"");
        hfTrueChecksCuts.push_back(*hfTrueChecksCutTemp);
        //---- Histos: Trigger info
        h_CBEsum[i] = hfTrigChecksCuts.at(i).makeTH1D(Form("CBEsum %s",cuttitle[i].c_str()),"CB Esum","",BinSettings(250,0.,2000.),Form("h_CBEsum%s",cutname[i].c_str()),true);
        //---- Histos: Tagger info
        h_TaggTimeuw[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger time  unweighted%s",cuttitle[i].c_str()),"Time","",BinSettings(800,-400,400),Form("h_TaggTimeuw%s",cutname[i].c_str()),true);
        h_TaggTimeww[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger time  with weights%s",cuttitle[i].c_str()),"Time","",BinSettings(800,-400,400),Form("h_TaggTimeww%s",cutname[i].c_str()),true);
        h_TaggcorTimeuw[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger corrected time  unweighted%s",cuttitle[i].c_str()),"Time","",BinSettings(800,-400,400),Form("h_TaggcorTimeuw%s",cutname[i].c_str()),true);
        h_TaggcorTimeww[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger corrected time  with weights%s",cuttitle[i].c_str()),"Time","",BinSettings(800,-400,400),Form("h_TaggcorTimeww%s",cutname[i].c_str()),true);
        h_TaggPhEnww[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger photon energy with weights %s",cuttitle[i].c_str()),energy_bins,"Photon energy","",Form("h_TaggPhEnww%s",cutname[i].c_str()),true);
        h_TaggPhEnuw[i] = hfTaggChecksCuts.at(i).makeTH1D(Form("Tagger photon energy unweighted %s",cuttitle[i].c_str()),energy_bins,"Photon energy","",Form("h_TaggPhEnuw%s",cutname[i].c_str()),true);
        //---- Histos: Reconstructed candidates
        h_IMeegReco[i] = hfCandChecksCuts.at(i).makeTH1D(Form("IMeeg %s",cuttitle[i].c_str()),"IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),Form("h_IMeegReco%s",cutname[i].c_str()),true);
        h_IMggReco[i] = hfCandChecksCuts.at(i).makeTH1D(Form("IMgg %s",cuttitle[i].c_str()),"IM(#gamma#gamma)","",BinSettings(250,0.,1000.),Form("h_IMggReco%s",cutname[i].c_str()),true);
        h_MMpReco[i] = hfCandChecksCuts.at(i).makeTH1D(Form("MMp %s",cuttitle[i].c_str()),"MM(p)","",BinSettings(250,400.,1400.),Form("h_MMpReco%s",cutname[i].c_str()),true);
        h_OpAngpphReco[i] = hfCandChecksCuts.at(i).makeTH1D(Form("OpAngpph %s",cuttitle[i].c_str()),"cos(#theta^{*}(p - #sum#gamma))","",BinSettings(100,-1.,1.),Form("h_OpAngpphReco%s",cutname[i].c_str()),true);
        h_PIDEvsT[i] = hfCandChecksCuts.at(i).makeTH2D(Form("PID E vs T %s",cuttitle[i].c_str()),"PID Energy","PID Time",BinSettings(200,0.,20),BinSettings(241,-60.5,60.5),Form("h_PIDEvsT%s",cutname[i].c_str()),true);
        h_TVetoEvsT[i] = hfCandChecksCuts.at(i).makeTH2D(Form("TAPSVeto E vs T %s",cuttitle[i].c_str()),"TAPSVeto Energy","TAPSVeto Time",BinSettings(200,0.,20),BinSettings(241,-60.5,60.5),Form("h_TVetoEvsT%s",cutname[i].c_str()),true);
        h_CBEvsT[i] = hfCandChecksCuts.at(i).makeTH2D(Form("CB E vs T %s",cuttitle[i].c_str()),"CB Energy","CB Time",BinSettings(200,0.,600),BinSettings(241,-60.5,60.5),Form("h_CBEvsT%s",cutname[i].c_str()),true);
        h_TAPSEvsT[i] = hfCandChecksCuts.at(i).makeTH2D(Form("TAPS E vs T %s",cuttitle[i].c_str()),"TAPS Energy","TAPS Time",BinSettings(200,0.,600),BinSettings(241,-60.5,60.5),Form("h_TAPSEvsT%s",cutname[i].c_str()),true);
        h_CBEvsNrCr[i] = hfCandChecksCuts.at(i).makeTH2D(Form("CB E vs NrCr %s",cuttitle[i].c_str()),"CB Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_CBEvsNrCr%s",cutname[i].c_str()),true);
        h_TAPSEvsNrCr[i] = hfCandChecksCuts.at(i).makeTH2D(Form("TAPS E vs NrCr %s",cuttitle[i].c_str()),"TAPS Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_TAPSEvsNrCr%s",cutname[i].c_str()),true);
        h_TimeCBvsPID[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Time CB vs PID %s",cuttitle[i].c_str()),"CB Time","PID Time",BinSettings(241,-60.5,60.5),BinSettings(241,-60.5,60.5),Form("h_TimeCBvsPID%s",cutname[i].c_str()),true);
        h_EnergyCBvsPID[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Energy CB vs PID %s",cuttitle[i].c_str()),"CB Energy","PID Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyCBvsPID%s",cutname[i].c_str()),true);
        h_TimeTAPSvsTVeto[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Time TAPS vs TAPSVeto %s",cuttitle[i].c_str()),"TAPS Time","TAPSVeto Time",BinSettings(241,-60.5,60.5),BinSettings(241,-60.5,60.5),Form("h_TimeTAPSvsTAPSVeto%s",cutname[i].c_str()),true);
        h_EnergyTAPSvsTVeto[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Energy TAPS vs TAPSVeto %s",cuttitle[i].c_str()),"TAPS Energy","TAPSVeto Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyTAPSvsTAPSVeto%s",cutname[i].c_str()),true);
        h_ThPhCBvsPID[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Theta-Phi CB vs PID %s",cuttitle[i].c_str()),"#theta_{CB-PID}","#phi_{CB-PID}",BinSettings(400,-100.,100.),BinSettings(400,-100.,100.),Form("h_ThPhCBvsPID%s",cutname[i].c_str()),true);
        h_ThPhTAPSvsTVeto[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Theta-Phi angle TAPS vs TAPSVeto %s",cuttitle[i].c_str()),"#theta_{TAPS-TAPSVeto}","#phi_{TAPS-TAPSVeto}",BinSettings(400,-100.,100.),BinSettings(400,-100.,100.),Form("h_ThPhTAPSvsTVeto%s",cutname[i].c_str()),true);
        h_ThetavsEnergy[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Theta vs energy %s",cuttitle[i].c_str()),"#theta [deg]","Cluster energy",BinSettings(180,-0.5,179.5),BinSettings(200,0.,600),Form("h_ThetavsEnergy%s",cutname[i].c_str()),true);
        for(int j=0; j<nrPartTypes; j++){
            h_EnergyCBvsPID_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("Energy CB vs PID %s for %s",cuttitle[i].c_str(),particletitle[j].c_str()),"CB Energy","PID Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyCBvsPID_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
            h_EnergyTAPSvsTVeto_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("Energy TAPS vs TAPSVeto %s for %s",cuttitle[i].c_str(),particletitle[j].c_str()),"TAPS Energy","TAPSVeto Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyTAPSvsTAPSVeto_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
            h_ThetavsEnergy_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("Theta vs energy %s for %s",cuttitle[i].c_str(),particletitle[j].c_str()),"#theta [deg]","Cluster energy",BinSettings(180,-0.5,179.5),BinSettings(200,0.,600),Form("h_ThetavsEnergy_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
            h_CBEvsNrCr_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("CB E vs NrCr %s for %s",cuttitle[i].c_str(),particletitle[j].c_str()),"CB Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_CBEvsNrCr_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
            h_TAPSEvsNrCr_RecMat[i][j] = hfMatchedTemp->makeTH2D(Form("TAPS E vs NrCr %s for %s",cuttitle[i].c_str(),particletitle[j].c_str()),"TAPS Energy","Nr crystals/cluster",BinSettings(200,0.,600),BinSettings(50,-0.5,49.5),Form("h_TAPSEvsNrCr_RecMat%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
        }
        h_EnergyCBvsPID_RecMat[i][nrPartTypes] = hfMatchedTemp->makeTH2D(Form("Energy CB vs PID %s for photon, only 1 PID hit",cuttitle[i].c_str()),"CB Energy","PID Energy",BinSettings(200,0.,600),BinSettings(200,0.,20),Form("h_EnergyCBvsPID_RecMat%sphoton1PIDhit",cutname[i].c_str()),true);
        h_OpAngpphReco_RecMat[i] = hfMatchedTemp->makeTH1D(Form("OpAngpph %s",cuttitle[i].c_str()),"cos(#theta^{*}(p - #sum#gamma))","",BinSettings(100,-1.,1.),Form("h_OpAngpphReco_RecMat%s",cutname[i].c_str()),true);
        h_NrRecCand[i] = hfCandChecksCuts.at(i).makeTH2D(Form("Nr reconstructed candidates %s",cuttitle[i].c_str()),"Nr cand in CB","Nr cand in TAPS",BinSettings(20,-0.5,19.5),BinSettings(20,-0.5,19.5),Form("h_NrRecCand%s",cutname[i].c_str()),true);
        //---- Histos: True MC
        h_ee_angle[i] = hfTrueChecksCuts.at(i).makeTH1D(Form("Cos opening angle e^{+}e^{-}, %s",cuttitle[i].c_str()),"Cos opening angle e^{+}e^{-}","",BinSettings(100,-1.,1.),Form("h_cosee_angle_%s",cutname[i].c_str()),true);
        h_gg_angle[i] = hfTrueChecksCuts.at(i).makeTH1D(Form("Cos opening angle #gamma#gamma, %s",cuttitle[i].c_str()),"Cos opening angle #gamma#gamma","",BinSettings(100,-1.,1.),Form("h_cosgg_angle_%s",cutname[i].c_str()),true);
        for(int j=0; j<nrPartTypes; j++){
            h_ThetavsEnergy_MCTrue[i][j] = hfTrueChecksCuts.at(i).makeTH2D(Form("Energy vs #theta for %s %s",particletitle[j].c_str(),cuttitle[i].c_str()),"#theta [deg]","Energy",BinSettings(180,-0.5,179.5),BinSettings(200,0.,600.),Form("h_ThetavsEnergy_MCTrue%s%s",cutname[i].c_str(),particlename[j].c_str()),true);
        }
        //---- delete the temporary histogramfactories
        delete hfTrigChecksCutTemp; delete hfTaggChecksCutTemp; delete hfCandChecksCutTemp; delete hfMatchedTemp; delete hfTrueChecksCutTemp;
    }
    //--- Overview
    h_AnalysisStat = hfOverview->makeTH1D("Analysis statistics","","",BinSettings(3+nrEvSel,-0.5,(3+nrEvSel)-0.5),"h_AnalysisStat",true);
    h_AnalysisStat_RecMat = hfOverview->makeTH2D("Analysis statistics","","",BinSettings(3+nrEvSel,-0.5,(3+nrEvSel)-0.5),BinSettings(nrPartTypes,-0.5,nrPartTypes-0.5),"h_AnalysisStat_RecMat",true);
    h_AnalysisStat->GetXaxis()->SetBinLabel(2,"Event"); h_AnalysisStat->GetXaxis()->SetBinLabel(3,"Valid tagghit");
    h_AnalysisStat_RecMat->GetXaxis()->SetBinLabel(1,"Event"); h_AnalysisStat_RecMat->GetXaxis()->SetBinLabel(2,"RecMat");
    h_AnalysisStat_RecMat->GetXaxis()->SetBinLabel(3,"Valid tagghit");
    for(int i=1; i<nrEvSel; i++){
        h_AnalysisStat->GetXaxis()->SetBinLabel(i+3,cutname[i].c_str());
        h_AnalysisStat_RecMat->GetXaxis()->SetBinLabel(i+3,cutname[i].c_str());
    }
    for(int i=0; i<nrPartTypes; i++)
        h_AnalysisStat_RecMat->GetYaxis()->SetBinLabel(i+1,particletitle[i].c_str());

    //--- KinFit
    string kfname[] = {"KF1p3ph","KF1p3phnp","KFnp3ph","KF1p2ph","KF1p2phnp","KFnp2ph"};
    string kftitle[] = {"KF 1p3ph","KF 1p3ph np","KF np3ph","KF 1p2ph","KF 1p2ph np","KF np2ph"};
    string fitvarnameCB[] = {"invEk","theta","phi","R"};
    string fitvarnameTA[] = {"invEk","Rxy","phi","L"};
    string fitvartitleCB[] = {"1/Ek","#theta","#phi","R"};
    string fitvartitleTA[] = {"1/Ek","R_{xy}","#phi","L"};
    vector<HistogramFactory> hfKinFits;
    for(int i=0; i<nrKF; i++){
        //---- Create subfolders for each KinFit
        auto hfKinFitTemp = new HistogramFactory(kfname[i],*hfKinFit,"");
        h_Stat[i] = hfKinFitTemp->makeTH1D(Form("Statistics, %s",kftitle[i].c_str()),"","",BinSettings(10,-0.5,9.5),Form("h_Stat_%s",kfname[i].c_str()),true);
        h_Stat[i]->GetXaxis()->SetBinLabel(1,"Before fit"); h_Stat[i]->GetXaxis()->SetBinLabel(2,"Success"); h_Stat[i]->GetXaxis()->SetBinLabel(3,">P");
        h_Stat[i]->GetXaxis()->SetBinLabel(4,"Succesful fit");
        h_Status[i] = hfKinFitTemp->makeTH1D(Form("KinFit status, %s",kftitle[i].c_str()),"","",BinSettings(7,-0.5,6.5),Form("h_Status_%s",kfname[i].c_str()),true);
        h_Status[i]->GetXaxis()->SetBinLabel(1,"Success"); h_Status[i]->GetXaxis()->SetBinLabel(2,"NoConvergence");  h_Status[i]->GetXaxis()->SetBinLabel(3,"TooManyIterations");
        h_Status[i]->GetXaxis()->SetBinLabel(4,"UnphysicalValues");  h_Status[i]->GetXaxis()->SetBinLabel(5,"NegativeDoF");  h_Status[i]->GetXaxis()->SetBinLabel(6,"OutOfMemory");
        h_Status[i]->GetXaxis()->SetBinLabel(7,"Unknown");
        h_Chi2[i] = hfKinFitTemp->makeTH1D(Form("#chi^{2}, %s",kftitle[i].c_str()),"","",BinSettings(100,-0.5,99.5),Form("h_Chi2_%s",kfname[i].c_str()),true);
        h_NDF[i] = hfKinFitTemp->makeTH1D(Form("NDF, %s",kftitle[i].c_str()),"","",BinSettings(100,-0.5,99.5),Form("h_NDF_%s",kfname[i].c_str()),true);
        h_EP[i] = hfKinFitTemp->makeTH2D(Form("E vs P, %s",kftitle[i].c_str()),"E_{initial}^{fit} - E_{final}^{fit}","P_{initial}^{fit} - P_{final}^{fit}",BinSettings(101,-200,200),BinSettings(101,-200,200),Form("h_EP_%s",cutname[i].c_str()),true);
        h_NIter[i] = hfKinFitTemp->makeTH2D(Form("Nr iterations per status, %s",kftitle[i].c_str()),"Nr iteratations","",BinSettings(100,0,100),BinSettings(7,-0.5,6.5),Form("h_NIter_%s",cutname[i].c_str()),true);
        h_NIter[i]->GetYaxis()->SetBinLabel(1,"Success"); h_NIter[i]->GetYaxis()->SetBinLabel(2,"NoConvergence");  h_NIter[i]->GetYaxis()->SetBinLabel(3,"TooManyIterations");
        h_NIter[i]->GetYaxis()->SetBinLabel(4,"UnphysicalValues");  h_NIter[i]->GetYaxis()->SetBinLabel(5,"NegativeDoF");  h_NIter[i]->GetYaxis()->SetBinLabel(6,"OutOfMemory");
        h_NIter[i]->GetYaxis()->SetBinLabel(7,"Unknown");
        h_NFuncCall[i] = hfKinFitTemp->makeTH2D(Form("Nr function calls per status, %s",kftitle[i].c_str()),"Nr function calls","",BinSettings(200,0,2000),BinSettings(7,-0.5,6.5),Form("h_NFuncCall_%s",cutname[i].c_str()),true);
        h_NFuncCall[i]->GetYaxis()->SetBinLabel(1,"Success"); h_NFuncCall[i]->GetYaxis()->SetBinLabel(2,"NoConvergence");  h_NFuncCall[i]->GetYaxis()->SetBinLabel(3,"TooManyIterations");
        h_NFuncCall[i]->GetYaxis()->SetBinLabel(4,"UnphysicalValues");  h_NFuncCall[i]->GetYaxis()->SetBinLabel(5,"NegativeDoF");  h_NFuncCall[i]->GetYaxis()->SetBinLabel(6,"OutOfMemory");
        h_NFuncCall[i]->GetYaxis()->SetBinLabel(7,"Unknown");
        for(int j=0; j<nrPcut; j++){
            //-- Create subfolders for each Pcut
            auto hfPcutstemp = new HistogramFactory(Form("Pcut%04d",(int)(pcuts[j]*1000)), *hfKinFitTemp);
            auto hfPullsPcutstempCB = new HistogramFactory("PullsCB", *hfPcutstemp);
            auto hfPullsPcutstempTA = new HistogramFactory("PullsTAPS", *hfPcutstemp);
            h_Prob[i][j] = hfPcutstemp->makeTH1D(Form("Probability, %s, pcut%04d",kftitle[i].c_str(),(int)(pcuts[j]*1000)),"P(#chi^{2})","",BinSettings(1000,0.,1.),Form("h_Prob_%s_pcut%04d",kfname[i].c_str(),(int)(pcuts[j]*1000)),true);
            h_Zv[i][j] = hfPcutstemp->makeTH1D(Form("Z-vertex, %s, pcut%04d",kftitle[i].c_str(),(int)(pcuts[j]*1000)),"Z vertex","",BinSettings(100,-15.,15.),Form("h_Zv_%s_pcut%04d",kfname[i].c_str(),(int)(pcuts[j]*1000)),true);
            h_IM3g[i][j] = hfPcutstemp->makeTH1D(Form("IM(e^{+}e^{-}#gamma), %s, pcut%04d",kftitle[i].c_str(),(int)(pcuts[j]*1000)),"IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),Form("h_IM3g_%s_pcut%04d",kfname[i].c_str(),(int)(pcuts[j]*1000)),true);
            h_IM2g[i][j] = hfPcutstemp->makeTH1D(Form("IM(#gamma#gamma), %s, pcut%04d",kftitle[i].c_str(),(int)(pcuts[j]*1000)),"IM(#gamma#gamma)","",BinSettings(250,0.,1000.),Form("h_IM2g_%s_pcut%04d",kfname[i].c_str(),(int)(pcuts[j]*1000)),true);
            h_EbeamPulls[i][j] = hfPcutstemp->makeTH1D(Form("Ebeam pull, %s, pcut%04d",kftitle[i].c_str(),(int)(pcuts[j]*1000)),"Pull","",BinSettings(200,-10.,10.),Form("h_EbeamPulls_%s_pcut%04d",kfname[i].c_str(),(int)(pcuts[j]*1000)),true);
            h_ZvertPulls[i][j] = hfPcutstemp->makeTH1D(Form("Z-vertex pull, %s, pcut%04d",kftitle[i].c_str(),(int)(pcuts[j]*1000)),"Pull","",BinSettings(200,-10.,10.),Form("h_ZvertPulls_%s_pcut%04d",kfname[i].c_str(),(int)(pcuts[j]*1000)),true);
            for(int k=0; k<nrPartTypes; k++){
                for(int l=0; l<nrKFVars; l++){
                    h_PartPulls_CB[i][j][k][l] = hfPullsPcutstempCB->makeTH1D(Form("%s pull for %s, %s, pcut%04d",fitvartitleCB[l].c_str(),particletitle[k].c_str(),kftitle[i].c_str(),(int)(pcuts[j]*1000)),Form("%s pull",fitvartitleCB[l].c_str()),"",BinSettings(200,-10.,10.),Form("h_PartPulls_%s_%s_%s_pcut%04d",fitvarnameCB[l].c_str(),particlename[k].c_str(),kfname[i].c_str(),(int)(pcuts[j]*1000)),true);
                    h_PartPulls_TAPS[i][j][k][l] = hfPullsPcutstempTA->makeTH1D(Form("%s pull for %s, %s, pcut%04d",fitvartitleTA[l].c_str(),particletitle[k].c_str(),kftitle[i].c_str(),(int)(pcuts[j]*1000)),Form("%s pull",fitvartitleTA[l].c_str()),"",BinSettings(200,-10.,10.),Form("h_PartPulls_%s_%s_%s_pcut%04d",fitvarnameTA[l].c_str(),particlename[k].c_str(),kfname[i].c_str(),(int)(pcuts[j]*1000)),true);
                }
            }
            delete hfPcutstemp; delete hfPullsPcutstempCB; delete hfPullsPcutstempTA;
        }
        delete hfKinFitTemp;
    }
}

void Pi0Dalitz::DoTrueMCStuff(const int cut, const std::vector<bool> &WhichMC, const std::vector<TParticlePtr> &trueparts, const double &tw)
{
    TParticlePtr pep; TParticlePtr pem;
    for(auto& truepart: trueparts){
        if(truepart->Type() == ParticleTypeDatabase::Proton) h_ThetavsEnergy_MCTrue[cut][en_p]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
        if(truepart->Type() == ParticleTypeDatabase::ePlus){
            h_ThetavsEnergy_MCTrue[cut][en_ep]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
            pep = truepart;
        }
        if(truepart->Type() == ParticleTypeDatabase::eMinus){
            h_ThetavsEnergy_MCTrue[cut][en_em]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
            pem = truepart;
        }
        if(truepart->Type() == ParticleTypeDatabase::Photon) h_ThetavsEnergy_MCTrue[cut][en_g]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
        //-- Additionally, if the MC channel is not eeg or gg, then fill neutrons-> protons, pi+->e+ and pi- -> e-
        if(WhichMC.at(2)){
            if(truepart->Type() == ParticleTypeDatabase::Neutron) h_ThetavsEnergy_MCTrue[cut][en_p]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
            if(truepart->Type() == ParticleTypeDatabase::PiPlus) h_ThetavsEnergy_MCTrue[cut][en_ep]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
            if(truepart->Type() == ParticleTypeDatabase::PiMinus) h_ThetavsEnergy_MCTrue[cut][en_em]->Fill(truepart->Theta()*radtodeg,truepart->Ek(),tw);
        }
    }
    //-- pi0->eeg
    if(WhichMC.at(0)){
        if(cut == 0){
            h_IMeegTrue->Fill((*trueparts.at(1)+*trueparts.at(2)+*trueparts.at(3)).M(),tw);
            h_IMeeTrue->Fill((*pep+*pem).M(),tw);
        }
        h_ee_angle[cut]->Fill(cos(TParticle::CalcAngle(pep,pem)),tw);
    }
    //-- pi0->gg
    if(WhichMC.at(1)){
        if(cut == 0) h_IMggTrue->Fill((*trueparts.at(1)+*trueparts.at(2)).M(),tw);
        h_gg_angle[cut]->Fill(cos(TParticle::CalcAngle(trueparts.at(1),trueparts.at(2))),tw);
    }
}

void Pi0Dalitz::DoCandSelStuff(const TCandidateList &recocands, TCandidatePtrList &selcands, std::vector<TParticlePtr> &selprot, std::vector<TParticleList> &selphots, std::vector<bool>& cases)
{
    TCandidatePtrList selectedCBcha, selectedCBneu, selectedTAPScha, selectedTAPSneu;
    for(const auto& currcand : recocands.get_iter()){
        bool inCB = false;
        if(currcand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB) inCB = true;
        //-- Before selection
        PlotCandSelStuff(currcand, 0);
        //-- Possibly tighter time cuts
        PlotCandSelStuff(currcand, 1);
        //-- CB candidates
        if(inCB){
            //-- Select neutral or charged particles and apply the energy correction
            if(currcand->VetoEnergy <= pidecut_chaneu){
                TCandidatePtr newcand = make_shared<TCandidate>(ParticleECorr(currcand,true,inCB));
                PlotCandSelStuff(newcand, 2);
                selcands.push_back(newcand);
                selectedCBneu.push_back(newcand);
            }
            else {
                TCandidatePtr newcand = make_shared<TCandidate>(ParticleECorr(currcand,false,inCB));
                PlotCandSelStuff(newcand, 3);
                selcands.push_back(newcand);
                selectedCBcha.push_back(newcand);
            }
        }
        else {
            //-- TAPS candidates, select neutral or charged particles and apply the energy correction (for neutral at least)
            if(currcand->VetoEnergy <= vetoecut_chaneu){
                TCandidatePtr newcand = make_shared<TCandidate>(ParticleECorr(currcand,true,inCB));
                PlotCandSelStuff(newcand, 4);
                selcands.push_back(newcand);
                selectedTAPSneu.push_back(newcand);
            }
            else {
                TCandidatePtr newcand = make_shared<TCandidate>(ParticleECorr(currcand,false,inCB));
                PlotCandSelStuff(newcand, 5);
                selcands.push_back(newcand);
                selectedTAPScha.push_back(newcand);
            }
        }
    }
    //-- Check that we only have the desired amount of selected candidates
    //   For Dalitz decay: 1 neutral CB, 2 charged CB and 0 or 1 proton candidate in TAPS or CB
    //   For 2g decay: 2 neutral CB and 0 or 1 proton candidate in TAPS or CB
    //   Additionally, sort the selected candidates into possible proton-photons combinations in preparation for the kinfit by applying certain
    //   "photon" and proton related requirements, currently only PID vs CB E-cut
    bool HasProtonCand = false;
    bool IsDalDecay = false;
    //-- Dalitz decay scenarios
    if(selectedTAPSneu.size() == 0 && ((selectedCBneu.size() == 1 && selectedCBcha.size() == 2 && selectedTAPScha.size() < 2) ||  (selectedCBneu.size() == 1 && selectedCBcha.size() == 3 && selectedTAPScha.size() == 0))){
        for(auto& cand : selcands){
            PlotCandSelStuff(cand,6);
        }
        //-- In the simple scenario of having only two charged candidates in CB, these are considered to be leptons (if they pass the E-cut)
        if(selectedCBcha.size() == 2){
            TParticleList tmpparts;
            tmpparts.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,selectedCBneu.at(0)));
            for(const auto& sellep : selectedCBcha){
                if(sellep->VetoEnergy < pidecut_lepcan)
                    tmpparts.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,sellep));
            }
            //-- If all lepton-candidates passed the checks, keep them as photon-candidates and check for a proton in TAPS
            if(tmpparts.size() == 3){
                //-- If we additionally have one candidate in TAPS, only this candidate is considered as the proton, if it pass the E-cut
                //   if it doesn't pass the E-cut then this is not a viable Dalitz decay event
                if(selectedTAPScha.size() == 1){
                    if(selectedTAPScha.at(0)->VetoEnergy && (selectedTAPScha.at(0)->VetoEnergy > (vetotacut_protcan_sl*selectedTAPScha.at(0)->CaloEnergy + vetotacut_protcan_in))){
                        selprot.push_back(make_shared<TParticle>(ParticleTypeDatabase::Proton,selectedTAPScha.at(0)));
                        selphots.push_back(tmpparts);
                        HasProtonCand = true;
                    }
                }
                else {
                    selphots.push_back(tmpparts);
                }
            }
        }
        //-- In the scenario of three charged candidates in CB, either one can be considered as the proton
        //   if they lie above a simple cut in E CB vs PID and the remainders are consered as leptons if they lie below another
        //   simple cut in E CB vs PID
        if(selectedCBcha.size() == 3){
            for(int i=0; i<3; i++){
                if(selectedCBcha.at(i)->VetoEnergy > (pidcbcut_protcan_sl*selectedCBcha.at(i)->CaloEnergy + pidcbcut_protcan_in)){
                    TParticleList tmpparts;
                    tmpparts.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,selectedCBneu.at(0)));
                    for(int j=0; j<3; j++){
                        if(i==j) continue;
                        if(selectedCBcha.at(j)->VetoEnergy < pidecut_lepcan)
                            tmpparts.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,selectedCBcha.at(j)));
                    }
                    //-- If all lepton-candidates passed the checks, keep them as photon-candidates and also store the proton-candidate
                    if(tmpparts.size() == 3){
                        selprot.push_back(make_shared<TParticle>(ParticleTypeDatabase::Proton,selectedCBcha.at(i)));
                        selphots.push_back(tmpparts);
                        HasProtonCand = true;
                    }
                }
            }
        }
        if(selphots.size() > 0){
            IsDalDecay = true;
            for(const auto& vpph : selphots){
                for(const auto& pph : vpph){
                    PlotCandSelStuff(pph->Candidate,7);
                }
            }
            for(const auto& ppr : selprot){
                PlotCandSelStuff(ppr->Candidate,8);
            }
        }
    }
    //-- 2g decay scenarios
    bool Is2gDecay = false;
    if(selectedTAPSneu.size() == 0 && ((selectedCBneu.size() == 2 && selectedCBcha.size() == 0 && selectedTAPScha.size() < 2) ||  (selectedCBneu.size() == 2 && selectedCBcha.size() == 1 && selectedTAPScha.size() == 0))){
        for(auto& cand : selcands){
            PlotCandSelStuff(cand,9);
        }
        TParticleList tmpparts;
        tmpparts.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,selectedCBneu.at(0)));
        tmpparts.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,selectedCBneu.at(1)));
        //-- if there's a charged candidate in CB or TAPS, check that it passes the proton E-cut. If not, it's not a viable 2g decay event
        if(selectedCBcha.size() == 1){
            if(selectedCBcha.at(0)->VetoEnergy > (pidcbcut_protcan_sl*selectedCBcha.at(0)->CaloEnergy + pidcbcut_protcan_in)){
                selprot.push_back(make_shared<TParticle>(ParticleTypeDatabase::Proton,selectedCBcha.at(0)));
                selphots.push_back(tmpparts);
                HasProtonCand = true;
            }
        }
        else if(selectedTAPScha.size() == 1){
            if(selectedTAPScha.at(0)->VetoEnergy && (selectedTAPScha.at(0)->VetoEnergy > (vetotacut_protcan_sl*selectedTAPScha.at(0)->CaloEnergy + vetotacut_protcan_in))){
                selprot.push_back(make_shared<TParticle>(ParticleTypeDatabase::Proton,selectedTAPScha.at(0)));
                selphots.push_back(tmpparts);
                HasProtonCand = true;
            }
        }
        else
            selphots.push_back(tmpparts);

        if(selphots.size() > 0){
            Is2gDecay = true;
            for(const auto& vpph : selphots){
                for(const auto& pph : vpph){
                    PlotCandSelStuff(pph->Candidate,10);
                }
            }
            for(const auto& ppr : selprot){
                PlotCandSelStuff(ppr->Candidate,11);
            }
        }
    }
    //-- Lastly store the checks of the different cases
    cases.push_back(HasProtonCand);
    cases.push_back(IsDalDecay);
    cases.push_back(Is2gDecay);
}

void Pi0Dalitz::PlotCandSelStuff(const TCandidatePtr cand, const int sel)
{
    TClusterPtr caloclu = cand->FindCaloCluster();
    TClusterPtr vetoclu = nullptr; if(cand->VetoEnergy) vetoclu = cand->FindVetoCluster();

    //-- Fill the calorimeter related info
    if(caloclu->DetectorType == Detector_t::Type_t::CB){
        h_CBEvsT_cs[sel]->Fill(cand->CaloEnergy,caloclu->Time);
        h_CBEvsNrCr_cs[sel]->Fill(cand->CaloEnergy,cand->ClusterSize);
        if(cand->VetoEnergy){
            h_TimeCBvsPID_cs[sel]->Fill(caloclu->Time,vetoclu->Time);
            h_EnergyCBvsPID_cs[sel]->Fill(cand->CaloEnergy,cand->VetoEnergy);
        }
    }
    else {
        h_TAPSEvsT_cs[sel]->Fill(cand->CaloEnergy,caloclu->Time);
        h_TAPSEvsNrCr_cs[sel]->Fill(cand->CaloEnergy,cand->ClusterSize);
        if(cand->VetoEnergy){
            h_TimeTAPSvsTVeto_cs[sel]->Fill(caloclu->Time,vetoclu->Time);
            h_EnergyTAPSvsTVeto_cs[sel]->Fill(cand->CaloEnergy,cand->VetoEnergy);
        }
    }
    // -- Fill the veto related info
    if(cand->VetoEnergy && vetoclu->DetectorType == Detector_t::Type_t::PID){
        h_PIDEvsT_cs[sel]->Fill(vetoclu->Energy,vetoclu->Time);
    }
    if(cand->VetoEnergy && vetoclu->DetectorType == Detector_t::Type_t::TAPSVeto){
        h_TVetoEvsT_cs[sel]->Fill(vetoclu->Energy,vetoclu->Time);
    }
}

TCandidate Pi0Dalitz::ParticleECorr(const TCandidatePtr cand, const bool &photon, const bool &incb)
{
    double origCluE = cand->CaloEnergy;
    double CluECorr = 0;
    if(doecorr){
        if(photon){
            if(incb)
                CluECorr = photCBcorr.GetECorr(origCluE);
            else
                CluECorr = photTAPScorr.GetECorr(origCluE);
        }
        else {
            if(incb)
                CluECorr = lepCBcorr.GetECorr(origCluE);
        }
    }
    double newCluE = origCluE + CluECorr;
    if(newCluE<0) newCluE=0;
    TClusterList cl;
    for_each(cand->Clusters.begin(), cand->Clusters.end(),[&cl](const TCluster& c){
        cl.emplace_back(c.Position,
                        c.Energy,
                        c.Time,
                        c.DetectorType,
                        c.CentralElement,
                        c.Hits);
                });
    TCandidate newcand = TCandidate(cand->Detector,
                         newCluE,
                         cand->Theta,
                         cand->Phi,
                         cand->Time,
                         cand->ClusterSize,
                         cand->VetoEnergy,
                         cand->TrackerEnergy,
                         cand->VetoEnergy > 0 ? TClusterList{std::prev(cl.end(), 2), std::prev(cl.end(), 1)}
                                             : TClusterList{std::prev(cl.end())}
                                               );
    return newcand;
}

void Pi0Dalitz::DoMatchTrueRecoStuff(const TParticleList &allmcpart, const std::vector<TParticlePtr> &trueparts, const TCandidatePtrList& recocands, std::vector<TParticlePtr> &matchrecopart)
{
    const auto matched  = utils::match1to1(allmcpart, recocands,
                                           [] (const TParticlePtr& p1, const TCandidatePtr& p2) {return p1->Angle(*p2);},
                                           {0.0, std_ext::degree_to_radian(15.0)});

    h_RecoTrueMatch->Fill(0);
    int nrgamma = 0; double Ekrecg=0; double Ektrueg=0;
    for(auto& truepart: trueparts){
        TCandidatePtr match = utils::FindMatched(matched,truepart);
        if(match){
            if(truepart->Type() == ParticleTypeDatabase::Proton){
                matchrecopart.push_back(make_shared<TParticle>(ParticleTypeDatabase::Proton,match));
                h_RecoTrueMatch->Fill(en_p+1);
                h_RecoTrueAngle->Fill(en_p+1,TParticle::CalcAngle(matchrecopart.back(),truepart)*radtodeg);
                if(matchrecopart.back()->Ek()>0) h_EktrueEkrec[en_p]->Fill(truepart->Ek()/matchrecopart.back()->Ek());
                h_AnalysisStat_RecMat->Fill(1.,en_p);
            }
            if(truepart->Type() == ParticleTypeDatabase::ePlus){
                matchrecopart.push_back(make_shared<TParticle>(ParticleTypeDatabase::ePlus,match));
                h_RecoTrueMatch->Fill(en_ep+1);
                h_RecoTrueAngle->Fill(en_ep+1,TParticle::CalcAngle(matchrecopart.back(),truepart)*radtodeg);
                if(matchrecopart.back()->Ek()>0) h_EktrueEkrec[en_ep]->Fill(truepart->Ek()/matchrecopart.back()->Ek());
                h_AnalysisStat_RecMat->Fill(1.,en_ep);
            }
            if(truepart->Type() == ParticleTypeDatabase::eMinus){
                matchrecopart.push_back(make_shared<TParticle>(ParticleTypeDatabase::eMinus,match));
                h_RecoTrueMatch->Fill(en_em+1);
                h_RecoTrueAngle->Fill(en_em+1,TParticle::CalcAngle(matchrecopart.back(),truepart)*radtodeg);
                if(matchrecopart.back()->Ek()>0) h_EktrueEkrec[en_em]->Fill(truepart->Ek()/matchrecopart.back()->Ek());
                h_AnalysisStat_RecMat->Fill(1.,en_em);
            }
            if(truepart->Type() == ParticleTypeDatabase::Photon){
                matchrecopart.push_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,match));
                h_RecoTrueMatch->Fill(en_g+1);
                h_RecoTrueAngle->Fill(en_g+1,TParticle::CalcAngle(matchrecopart.back(),truepart)*radtodeg);
                if(matchrecopart.back()->Ek()>0) h_EktrueEkrec[en_g]->Fill(truepart->Ek()/matchrecopart.back()->Ek());
                h_AnalysisStat_RecMat->Fill(1.,en_g);
                nrgamma++;
                Ekrecg+=matchrecopart.back()->Ek();
                Ektrueg+=truepart->Ek();
            }
        }
    }
    if(nrgamma==2 && Ekrecg>0)
        h_EktrueEkrec_gg->Fill(Ektrueg/Ekrecg);
}

void Pi0Dalitz::DoTaggerStuff(const int cut, const TLorentzVector &g, const double &time, const double &cortime, const double &tw)
{
    h_TaggTimeuw[cut]->Fill(time);
    h_TaggTimeww[cut]->Fill(time,tw);
    h_TaggcorTimeuw[cut]->Fill(cortime);
    h_TaggcorTimeww[cut]->Fill(cortime,tw);
    h_TaggPhEnww[cut]->Fill(g.E(),tw);
    h_TaggPhEnuw[cut]->Fill(g.E());

}

void Pi0Dalitz::DoTriggerStuff(const int cut, const double &tw)
{
    h_CBEsum[cut]->Fill(triggersimu.GetCBEnergySum(),tw);
}

void Pi0Dalitz::DoRecoCandStuff(const int cut, const TCandidatePtrList& recocands, const std::vector<TParticlePtr>& selprot, const std::vector<TParticleList>& selphots, const std::vector<TParticlePtr> &recmatparts, const std::vector<bool> &WhichMC, const TLorentzVector &ig, const double &tw)
{
    h_AnalysisStat->Fill(2+cut,tw);

    //-- Fill the reconstructed candidate information
    std::map<int, int> PIDElFreq;
    std::map<int, int> TVetoElFreq;
    int nrcandCB=0; int nrcandTAPS=0;
    for(const auto& currcand : recocands) {
        TClusterPtr caloclu = currcand->FindCaloCluster();
        TClusterPtr vetoclu = nullptr; if(currcand->VetoEnergy) vetoclu = currcand->FindVetoCluster();
        //--- PID/Veto timings, only once per veto hit, i.e. exclude duplicates
        //    also check which pid and tapsveto clusters are used multiple times
        //    two maps are created which keeps track of how often a veto element is used
        if(currcand->VetoEnergy){
            int elnr = vetoclu->CentralElement;
            if(vetoclu->DetectorType == Detector_t::Type_t::PID){
                auto result = PIDElFreq.insert(std::pair<int, int>(elnr, 1));
                if (result.second == false)
                    result.first->second++;
                else
                    h_PIDEvsT[cut]->Fill(vetoclu->Energy,vetoclu->Time,tw);
            }
            else{
                auto result = TVetoElFreq.insert(std::pair<int, int>(elnr, 1));
                if (result.second == false)
                    result.first->second++;
                else
                    h_TVetoEvsT[cut]->Fill(vetoclu->Energy,vetoclu->Time,tw);
            }
        }
        //--- CB related info
        if(caloclu->DetectorType == Detector_t::Type_t::CB){
            nrcandCB++;
            h_CBEvsT[cut]->Fill(currcand->CaloEnergy,caloclu->Time,tw);
            if(currcand->VetoEnergy && vetoclu->DetectorType == Detector_t::Type_t::PID){
                h_EnergyCBvsPID[cut]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
                h_TimeCBvsPID[cut]->Fill(currcand->Time,vetoclu->Time,tw);
                h_ThPhCBvsPID[cut]->Fill((caloclu->Position.Theta()-vetoclu->Position.Theta())*radtodeg,(caloclu->Position.Phi()-vetoclu->Position.Phi())*radtodeg,tw);
            }
            h_CBEvsNrCr[cut]->Fill(currcand->CaloEnergy,currcand->ClusterSize,tw);
        }
        //--- TAPS related info
        else{
            nrcandTAPS++;
            h_TAPSEvsT[cut]->Fill(currcand->CaloEnergy,currcand->Time,tw);
            //---- Here I'm ignoring the possibility (is it even possible, I don't think so) that a TAPS cluster can be paired with a PID cluster
            if(currcand->VetoEnergy && vetoclu->DetectorType == Detector_t::Type_t::TAPSVeto){
                h_EnergyTAPSvsTVeto[cut]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
                h_TimeTAPSvsTVeto[cut]->Fill(currcand->Time,vetoclu->Time,tw);
                h_ThPhTAPSvsTVeto[cut]->Fill((caloclu->Position.Theta()-vetoclu->Position.Theta())*radtodeg,(caloclu->Position.Phi()-vetoclu->Position.Phi())*radtodeg,tw);
            }
            h_TAPSEvsNrCr[cut]->Fill(currcand->CaloEnergy,currcand->ClusterSize,tw);
        }
        h_ThetavsEnergy[cut]->Fill(currcand->Theta*radtodeg,currcand->CaloEnergy,tw);
    }
    h_NrRecCand[cut]->Fill(nrcandCB,nrcandTAPS,tw);

    //-- How often are each PID/Veto element used per event
    if(cut==0){
        for(auto elem : PIDElFreq)
            h_PIDMultUsed->Fill(elem.second,elem.first,tw);
        for(auto elem : TVetoElFreq)
            h_VetoMultUsed->Fill(elem.second,elem.first,tw);
    }

    //-- Fill the info from the proton-photon particle combinations
    TLorentzVector InitialProtonVec = LorentzVec(vec3(0,0,0),ParticleTypeDatabase::Proton.Mass());
    TLorentzVector InitialStateVec = InitialProtonVec + ig;
    TVector3 InStBoost = -InitialStateVec.BoostVector();
    if(selprot.size()>0){
        //--- Fill combinatorially the opening angle between the proton and the photons in CM
        //    for this, at least one proton and one photon is needed
        for(int i=0; i<(int)selprot.size(); i++){
            TLorentzVector TLproton = *selprot.at(i);
            TLorentzVector TLallg;
            for(const auto& ph : selphots.at(i))
                TLallg += *ph;
            //---- Boost the proton and the summed photons vector to CM and check their opening angle
            TLproton.Boost(InStBoost);
            TLallg.Boost(InStBoost);
            h_OpAngpphReco[cut]->Fill(cos(TLproton.Angle(TLallg.Vect())),tw);
        }
    }
    //-- Using the photon-to-fit candidates, fill combinatorially the MMp, IMgg and IMggg histograms
    for(const auto& phots : selphots){
        TLorentzVector TLallg;
        for(const auto& ph : phots)
            TLallg += *ph;
        h_MMpReco[cut]->Fill((InitialProtonVec + ig - TLallg).M(),tw);
        if(phots.size() == 2)
            h_IMggReco[cut]->Fill(TLallg.M(),tw);
        if(phots.size() == 3)
            h_IMeegReco[cut]->Fill(TLallg.M(),tw);
    }

    //-- Fill histograms for the candidates matched with true MC
    TLorentzVector matvec_p;
    vector<TLorentzVector> matvec_pi0;
    int pt;
    for(auto& recmatpart: recmatparts){
        //--- What kind of particle is it
        if(recmatpart->Type() == ParticleTypeDatabase::Proton){pt = en_p; matvec_p = *recmatpart;}
        else if(recmatpart->Type() == ParticleTypeDatabase::ePlus){ pt = en_ep; matvec_pi0.push_back(*recmatpart);}
        else if(recmatpart->Type() == ParticleTypeDatabase::eMinus){ pt = en_em; matvec_pi0.push_back(*recmatpart);}
        else if (recmatpart->Type() == ParticleTypeDatabase::Photon){ pt = en_g; matvec_pi0.push_back(*recmatpart);}
        else continue;
        h_AnalysisStat_RecMat->Fill(2+cut,pt,tw);
        TCandidatePtr currcand = recmatpart->Candidate;
        h_ThetavsEnergy_RecMat[cut][pt]->Fill(radtodeg*currcand->Theta,currcand->CaloEnergy,tw);
        //--- CB related info
        if(currcand->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            h_CBEvsNrCr_RecMat[cut][pt]->Fill(currcand->CaloEnergy,currcand->ClusterSize,tw);
            if(currcand->VetoEnergy && currcand->FindVetoCluster()->DetectorType == Detector_t::Type_t::PID){
                h_EnergyCBvsPID_RecMat[cut][pt]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
                if(PIDElFreq.at(currcand->FindVetoCluster()->CentralElement) == 1 && pt == en_g)
                    h_EnergyCBvsPID_RecMat[cut][nrPartTypes]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
            }
        }
        //--- TAPS related info
        else{
            h_TAPSEvsNrCr_RecMat[cut][pt]->Fill(currcand->CaloEnergy,currcand->ClusterSize,tw);
            if(currcand->VetoEnergy && currcand->FindVetoCluster()->DetectorType == Detector_t::Type_t::TAPSVeto){
                h_EnergyTAPSvsTVeto_RecMat[cut][pt]->Fill(currcand->CaloEnergy,currcand->VetoEnergy,tw);
            }
        }
    }
    //--- Boost to CM frame and check the proton-pi0 opening angle
    if(((matvec_pi0.size()==3 && WhichMC.at(0)) || (matvec_pi0.size()==2 && WhichMC.at(1))) && matvec_p.E()){
        TLorentzVector summatvec_pi0;
        for(auto& matvec: matvec_pi0) summatvec_pi0 += matvec;
        summatvec_pi0.Boost(InStBoost);
        matvec_p.Boost(InStBoost);
        h_OpAngpphReco_RecMat[cut]->Fill(cos(matvec_p.Angle(summatvec_pi0.Vect())),tw);
    }
}

//-- This version requires a reconstructed proton
double Pi0Dalitz::DoKinFitStuffWProt(const int KFind, const std::vector<TParticlePtr>& selprot, const std::vector<TParticleList>& selphots, const TLorentzVector &ig, utils::KinFitter &fitobj, std::vector<TParticlePtr> &bestprobrec, std::vector<TParticlePtr> &bestprobfit, const double tw)
{
    TParticleList bestfitrecph, bestfitph;
    TParticlePtr bestfitrecp, bestfitp;
    std::vector<utils::Fitter::FitParticle> fitparticles;
    double fit_z_vert = -50;
    double fitbeamE = -50;
    double fitprob = -50;
    double pullbeamE = -50;
    double pullzvert = -50;
    double chi2 = -50;
    double ndf = -50;
    //-- loop over the proton-photon combinations
    for(int i=0; i<(int)selphots.size(); i++) {
        h_Stat[KFind]->Fill(0.,tw);
        TParticlePtr protontofit = selprot.at(i);
        //-- do the kinfit
        const auto& result = fitobj.DoFit(ig.E(), protontofit, selphots.at(i));
        //-- fill info related to the status of the fit
        h_Status[KFind]->Fill((int)result.Status,tw);
        h_NIter[KFind]->Fill(result.NIterations,(int)result.Status,tw);
        h_NFuncCall[KFind]->Fill(result.NFunctionCalls,(int)result.Status,tw);
        //-- if the fit did not converge or did not result in a better probability, try next combination
        if(result.Status != APLCON::Result_Status_t::Success)
            continue;
        h_Stat[KFind]->Fill(1.,tw);
        if(!std_ext::copy_if_greater(fitprob, result.Probability))
            continue;
        h_Stat[KFind]->Fill(2.,tw);

        //-- else, store the fit result and the reconstructed input
        bestfitrecp = protontofit;
        bestfitrecph = selphots.at(i);
        bestfitp = fitobj.GetFittedProton();
        bestfitph = fitobj.GetFittedPhotons();
        fitparticles = fitobj.GetFitParticles();
        fit_z_vert = fitobj.GetFittedZVertex();
        fitbeamE = fitobj.GetFittedBeamE();
        pullbeamE = fitobj.GetBeamEPull();
        pullzvert = fitobj.GetZVertexPull();
        chi2 = result.ChiSquare;
        ndf = result.NDoF;
    }

    if(fitprob==(-50)) return fitprob;
    h_Stat[KFind]->Fill(3., tw);
    //-- fill histograms
    TLorentzVector fitphsum; for(const auto& bfphot : bestfitph) fitphsum = fitphsum + *bfphot;
    TLorentzVector fitinitial = LorentzVec({0, 0, fitbeamE}, fitbeamE) + LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass());
    TLorentzVector IFdiff = fitinitial - *bestfitp - fitphsum;
    h_EP[KFind]->Fill(IFdiff.E(),IFdiff.P(),tw);
    h_Chi2[KFind]->Fill(chi2,tw);
    h_NDF[KFind]->Fill(ndf,tw);
    for(int i=0; i<nrPcut; i++){
        if(fitprob<pcuts[i]) break;
        h_Prob[KFind][i]->Fill(fitprob,tw);
        h_Zv[KFind][i]->Fill(fit_z_vert,tw);
        h_EbeamPulls[KFind][i]->Fill(pullbeamE,tw);
        h_ZvertPulls[KFind][i]->Fill(pullzvert,tw);
        if(bestfitph.size() == 3) h_IM3g[KFind][i]->Fill(fitphsum.M(),tw);
        if(bestfitph.size() == 2) h_IM2g[KFind][i]->Fill(fitphsum.M(),tw);
        //-- fill the pull histos depending on particle type and CB or TAPS, I can't distinguish between e+ and e- so it all goes in e-
        for(const auto &fp : fitparticles){
            bool inCB = false;
            if(fp.Particle->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB) inCB = true;
            int partype = -1;
            if(fp.Particle->Type() == ParticleTypeDatabase::Proton) partype = en_p;
            else if(fp.Particle->Type() == ParticleTypeDatabase::Photon){
                if((inCB && fp.Particle->Candidate->VetoEnergy <= pidecut_chaneu) || (!inCB && fp.Particle->Candidate->VetoEnergy <= vetoecut_chaneu))
                    partype = en_g;
                else
                    partype = en_em;
            }
            else {
                LOG(INFO)<<"In the fit particles was not identified as a proton or a photon";
            }
            for(int j=0; j<nrKFVars; j++){
                if(inCB)
                    h_PartPulls_CB[KFind][i][partype][j]->Fill(fp.GetPulls().at(j),tw);
                else
                    h_PartPulls_TAPS[KFind][i][partype][j]->Fill(fp.GetPulls().at(j),tw);
            }
        }
    }
    //-- fill the vectors with the rec and fit particles from the best fit
    for(const auto& ph: bestfitrecph)
        bestprobrec.push_back(ph);
    for(const auto& ph: bestfitph)
        bestprobfit.push_back(ph);
    bestprobrec.push_back(bestfitrecp);
    bestprobfit.push_back(bestfitp);

    return fitprob;
}

//-- This version doesn't require a reconstructed proton and in the current version of the analysis, there's only one list of reconstructed photons. i.e. no need to loop
double Pi0Dalitz::DoKinFitStuffNProt(const int KFind, const std::vector<TParticleList> &selphots, const TLorentzVector &ig, utils::NoProtonFitter& fitobj, TParticleList &fitphotons, LorentzVec &fitproton, const double tw)
{
    std::vector<utils::Fitter::FitParticle> fitparticles;
    double fit_z_vert = -50;
    double fitbeamE = -50;
    double fitprob = -50;
    double pullbeamE = -50;
    double pullzvert = -50;
    double chi2 = -50;
    double ndf = -50;
    h_Stat[KFind]->Fill(0.,tw);
    //-- check that there's indeed only one set of photons to fit
    if(selphots.size() != 1){
        LOG(INFO)<<"In the no-proton kinfit, there should only be one list of photons";
        return fitprob;
    }
    //-- do the kinfit
    const auto& result = fitobj.DoFit(ig.E(), selphots.at(0));
    //-- fill info related to the status of the fit
    h_Status[KFind]->Fill((int)result.Status,tw);
    h_NIter[KFind]->Fill(result.NIterations,(int)result.Status,tw);
    h_NFuncCall[KFind]->Fill(result.NFunctionCalls,(int)result.Status,tw);
    //-- check if the fit converged
    if(result.Status != APLCON::Result_Status_t::Success)
        return fitprob;
    h_Stat[KFind]->Fill(3.,tw);
    fitprob = result.Probability;
    //-- store the fit result
    fitproton = fitobj.GetFittedProton();
    fitphotons = fitobj.GetFittedPhotons();
    fitparticles = fitobj.GetFitParticles();
    fit_z_vert = fitobj.GetFittedZVertex();
    fitbeamE = fitobj.GetFittedBeamE();
    pullbeamE = fitobj.GetBeamEPull();
    pullzvert = fitobj.GetZVertexPull();
    chi2 = result.ChiSquare;
    ndf = result.NDoF;

    //-- fill histograms
    TLorentzVector fitphsum; for(const auto& bfphot : fitphotons) fitphsum = fitphsum + *bfphot;
    TLorentzVector fitinitial = LorentzVec({0, 0, fitbeamE}, fitbeamE) + LorentzVec({0,0,0}, ParticleTypeDatabase::Proton.Mass());
    TLorentzVector IFdiff = fitinitial - fitproton - fitphsum;
    h_EP[KFind]->Fill(IFdiff.E(),IFdiff.P(),tw);
    h_Chi2[KFind]->Fill(chi2,tw);
    h_NDF[KFind]->Fill(ndf,tw);
    for(int i=0; i<nrPcut; i++){
        if(fitprob<pcuts[i]) break;
        h_Prob[KFind][i]->Fill(fitprob,tw);
        h_Zv[KFind][i]->Fill(fit_z_vert,tw);
        h_EbeamPulls[KFind][i]->Fill(pullbeamE,tw);
        h_ZvertPulls[KFind][i]->Fill(pullzvert,tw);
        if(fitphotons.size() == 3) h_IM3g[KFind][i]->Fill(fitphsum.M(),tw);
        if(fitphotons.size() == 2) h_IM2g[KFind][i]->Fill(fitphsum.M(),tw);
        //-- fill the pull histos depending on particle type and CB or TAPS, I can't distinguish between e+ and e- so it all goes in e-
        for(const auto &fp : fitparticles){
            bool inCB = false;
            if(fp.Particle->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB) inCB = true;
            int partype = -1;
            if(fp.Particle->Type() == ParticleTypeDatabase::Proton) partype = en_p;
            else if(fp.Particle->Type() == ParticleTypeDatabase::Photon){
                if((inCB && fp.Particle->Candidate->VetoEnergy <= pidecut_chaneu) || (!inCB && fp.Particle->Candidate->VetoEnergy <= vetoecut_chaneu))
                    partype = en_g;
                else
                    partype = en_em;
            }
            else {
                LOG(INFO)<<"In the fit particles was not identified as a proton or a photon";
            }
            for(int j=0; j<nrKFVars; j++){
                if(inCB)
                    h_PartPulls_CB[KFind][i][partype][j]->Fill(fp.GetPulls().at(j),tw);
                else
                    h_PartPulls_TAPS[KFind][i][partype][j]->Fill(fp.GetPulls().at(j),tw);
            }
        }
    }

    return fitprob;
}

APLCON::Fit_Settings_t Pi0Dalitz::MakeFitSettings(unsigned max_iterations)
{
    APLCON::Fit_Settings_t settings;
    settings.MaxIterations = max_iterations;
    //settings.ConstraintAccuracy = 1.0e-3;
    //settings.Chi2Accuracy = 1.0e-2;
    //settings.DebugLevel = 0;
    return settings;
}

AUTO_REGISTER_PHYSICS(Pi0Dalitz)

