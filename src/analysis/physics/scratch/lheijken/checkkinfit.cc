#include "checkkinfit.h"

#include "expconfig/ExpConfig.h"
#include "utils/uncertainties/Interpolated.h"
#include "analysis/utils/uncertainties/FitterSergey.h"
#include "utils/ProtonPhotonCombs.h"
#include "utils/ParticleTools.h"

#include "base/Logger.h"

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
    fitter(nullptr, opts->Get<bool>("FitZVertex", true))
{
    fitter.SetZVertexSigma(3.0);

    CreateHistos();

    taps_detector = ExpConfig::Setup::GetDetector<expconfig::detector::TAPS>();
    TAPSZPos = taps_detector->GetZPosition();
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();
    CBRad = cb_detector->GetInnerRadius();
}

void scratch_lheijken_checkkinfit::ProcessEvent(const TEvent& event, manager_t&)
{
    if(!triggersimu.ProcessEvent(event))
        h_Steps->Fill("Triggersimu failed", 1.0);

    if(!triggersimu.HasTriggered())
        return;

    fitter.SetUncertaintyModel(fit_model);

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
    int nrfsparticles = 3;
    int nrphotonsinfit = 3;
    //-- Also at the moment, when analysing MC, let the decay type decide
    if(MCpi0eeg){ nrfsparticles = 3; nrphotonsinfit = 3; }
    if(MCpi0gg){ nrfsparticles = 2; nrphotonsinfit = 2; }
    if(MC2pi04g){ nrfsparticles = 4; nrphotonsinfit = 4; }
    if(nrcands < nrfsparticles) return;

    //-- Make all possible combinations of the candidates into groups containing one proton and the rest photons
    utils::ProtonPhotonCombs proton_photons(candidates);

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
    }

    //-- Tagger loop
    for(const TTaggerHit& taggerhit : data.TaggerHits) {

        h_Steps->Fill("Seen taggerhits",1.0);

        //-- Only bother with tagger hits which are inside the prompt or random windows
        promptrandom.SetTaggerTime(triggersimu.GetCorrectedTaggerTime(taggerhit));
        if(promptrandom.State() == PromptRandom::Case::Outside)
            continue;
        double taggw = promptrandom.FillWeight();

        //-- compare Ebeam True-Rec
        TParticlePtr truebeam = 0;
        for(auto& truepart: TruePart){
            if(truepart->Type() == ParticleTypeDatabase::BeamTarget){
                truebeam = truepart;
                h_TrueRec_Ebeam->Fill(taggerhit.PhotonEnergy,truebeam->Ek()-taggerhit.PhotonEnergy,taggw);
                break;
            }
        }

        //-- Do a kinematic fit for a fixed number of photons and 1 proton (only energy unmeasured)

        //-- create the 1 proton - X photons combinations, require the MM(Xg) to be close to proton mass
        auto filtered_combs = proton_photons()
                              .Observe([this] (const string& cut) { h_Steps->Fill(cut.c_str(), 1.0); }, "F ")
                              .FilterMM(taggerhit, ParticleTypeDatabase::Proton.GetWindow(350).Round());
        if(filtered_combs.empty()) {
            h_Steps->Fill("No combs left",1.0);
            continue;
        }
        TParticleList recphotonsbestfit; TParticleList fitphotonsbestfit;
        TParticlePtr recprotonbestfit = 0; TParticlePtr fitprotonbestfit = 0;
        double fitprob = std_ext::NaN;
        double fit_z_vert = -50;
        double fitbeamE = -50;
        //-- loop over the (filtered) proton combinations -- currently only if at least #cand >= # all f.s. particles
        if(nrcands > nrfsparticles){
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
                    if(!std_ext::copy_if_greater(fitprob, result.Probability))
                        continue;
                    //-- Else, store the fit result and the reconstructed input
                    recphotonsbestfit = photonstofit;
                    fitphotonsbestfit = fitter.GetFittedPhotons();
                    recprotonbestfit = protontofit;
                    fitprotonbestfit = fitter.GetFittedProton();
                    fit_z_vert = fitter.GetFittedZVertex();
                    fitbeamE = fitter.GetFittedBeamE();
                }
            }
            h_Probability->Fill(fitprob,taggw);
            if(fitprob > 0.001) {
                h_Steps->Fill("TreeFill",1.0);
                //-- Make all the comparisons between the protons and photons
                DoFitComparisons(fitprotonbestfit,recprotonbestfit,fitphotonsbestfit,recphotonsbestfit,TrueRecoMatchPart,TruePart,fit_z_vert,taggw);
                //-- Also for the beam and z-vertex
                if(truebeam){
                    h_TrueFit_Ebeam->Fill(truebeam->Ek(),truebeam->Ek() - fitbeamE,taggw);
                    h_TrueFit_zvert->Fill(true_z_vert,true_z_vert - fit_z_vert,taggw);
                }
                h_FitRec_Ebeam->Fill(taggerhit.PhotonEnergy, fitbeamE - taggerhit.PhotonEnergy,taggw);
                //-- Fill the IM histograms
                if(nrphotonsinfit == 3){
                    h_IMeeg_Fit->Fill((*fitphotonsbestfit.at(0) + *fitphotonsbestfit.at(1) + *fitphotonsbestfit.at(2)).M(),taggw);
                    h_IMeeg_Rec->Fill((*recphotonsbestfit.at(0) + *recphotonsbestfit.at(1) + *recphotonsbestfit.at(2)).M(),taggw);
                }
                if(nrphotonsinfit == 2){
                    h_IMgg_Fit->Fill((*fitphotonsbestfit.at(0) + *fitphotonsbestfit.at(1)).M(),taggw);
                    h_IMgg_Rec->Fill((*recphotonsbestfit.at(0) + *recphotonsbestfit.at(1)).M(),taggw);
                }
                if(nrphotonsinfit == 4){
                    for(int i=0; i<3; i++){
                        for(int j=(i+1);j<4; j++){
                            h_IMgg_Fit->Fill((*fitphotonsbestfit.at(i) + *fitphotonsbestfit.at(j)).M(),taggw);
                            h_IMgg_Rec->Fill((*recphotonsbestfit.at(i) + *recphotonsbestfit.at(j)).M(),taggw);
                        }
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
    auto hfTrueRecCB = new HistogramFactory("CB",*hfTrueRec,"");
    auto hfTrueRecTAPS = new HistogramFactory("TAPS",*hfTrueRec,"");
    auto hfFitRec = new HistogramFactory("FitRec",HistFac,"");
    auto hfFitRecCB = new HistogramFactory("CB",*hfFitRec,"");
    auto hfFitRecTAPS = new HistogramFactory("TAPS",*hfFitRec,"");
    auto hfTrueFit = new HistogramFactory("TrueFit",HistFac,"");
    auto hfTrueFitCB = new HistogramFactory("CB",*hfTrueFit,"");
    auto hfTrueFitTAPS = new HistogramFactory("TAPS",*hfTrueFit,"");

    //-- Histogram binnings and names
    string varname[] = {"Ek","theta","phi"};
    string vartitle[] = {"Ek","#theta","#phi"};
    string partname[] = {"p","ep","em","g"};
    string parttitle[] = {"p","e^{+}","e^{-}","#gamma"};
    const BinSettings EkBins = BinSettings(500,0,1000);
    const BinSettings ThCBBins = BinSettings(90,0,180);
    const BinSettings ThTABins = BinSettings(40,0,40);
    const BinSettings PhBins = BinSettings(200,-200,200);
    vector<BinSettings> VarBins;
    VarBins.push_back(BinSettings(101,-100,100));
    VarBins.push_back(BinSettings(101,-100,100));
    VarBins.push_back(BinSettings(101,-100,100));

    //-- The histograms
    h_Steps = hfOverview->makeTH1D("Steps","","",BinSettings(10),"h_Steps");
    h_Probability = hfOverview->makeTH1D("Fit probability","P(#chi^{2})","",BinSettings(1000,0.,1.),"h_Probability",true);
    //--- The kinematical variables for the final state particles
    for(int i=0; i<nrPartTypes; i++){
        for(int j=0; j<nrKinVars; j++){
            //-- TrueRec
            h_TrueRec_vsEk_CB[i][j] = hfTrueRecCB->makeTH2D(Form("TrueRec %s vs Ek^{rec} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"Ek [MeV]",Form("True - Rec %s",vartitle[j].c_str()),EkBins,VarBins.at(j),Form("hTrueRec%svsEkCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueRec_vsEk_TAPS[i][j] = hfTrueRecTAPS->makeTH2D(Form("TrueRec %s vs Ek^{rec} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"Ek [MeV]",Form("True - Rec %s",vartitle[j].c_str()),EkBins,VarBins.at(j),Form("hTrueRec%svsEkTAPS%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueRec_vsTh_CB[i][j] = hfTrueRecCB->makeTH2D(Form("TrueRec %s vs Th^{rec} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"#theta [deg]",Form("True - Rec %s",vartitle[j].c_str()),ThCBBins,VarBins.at(j),Form("hTrueRec%svsThCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueRec_vsTh_TAPS[i][j] = hfTrueRecTAPS->makeTH2D(Form("TrueRec %s vs Th^{rec} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"#theta [deg]",Form("True - Rec %s",vartitle[j].c_str()),ThTABins,VarBins.at(j),Form("hTrueRec%svsThTAPS%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueRec_vsPh_CB[i][j] = hfTrueRecCB->makeTH2D(Form("TrueRec %s vs Ph^{rec} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"#phi [deg]",Form("True - Rec %s",vartitle[j].c_str()),PhBins,VarBins.at(j),Form("hTrueRec%svsPhCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueRec_vsPh_TAPS[i][j] = hfTrueRecTAPS->makeTH2D(Form("TrueRec %s vs Ph^{rec} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"#phi [deg]",Form("True - Rec %s",vartitle[j].c_str()),PhBins,VarBins.at(j),Form("hTrueRec%svsPhTAPS%s",varname[j].c_str(),partname[i].c_str()),true);
            //-- FitRec
            h_FitRec_vsEk_CB[i][j] = hfFitRecCB->makeTH2D(Form("FitRec %s vs Ek^{rec} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"Ek [MeV]",Form("Fit - Rec %s",vartitle[j].c_str()),EkBins,VarBins.at(j),Form("hFitRec%svsEkCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_FitRec_vsEk_TAPS[i][j] = hfFitRecTAPS->makeTH2D(Form("FitRec %s vs Ek^{rec} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"Ek [MeV]",Form("Fit - Rec %s",vartitle[j].c_str()),EkBins,VarBins.at(j),Form("hFitRec%svsEkTAPS%s",varname[j].c_str(),partname[i].c_str()),true);
            h_FitRec_vsTh_CB[i][j] = hfFitRecCB->makeTH2D(Form("FitRec %s vs Th^{rec} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"#theta [deg]",Form("Fit - Rec %s",vartitle[j].c_str()),ThCBBins,VarBins.at(j),Form("hFitRec%svsThCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_FitRec_vsTh_TAPS[i][j] = hfFitRecTAPS->makeTH2D(Form("FitRec %s vs Th^{rec} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"#theta [deg]",Form("Fit - Rec %s",vartitle[j].c_str()),ThTABins,VarBins.at(j),Form("hFitRec%svsThTAPS%s",varname[j].c_str(),partname[i].c_str()),true);
            h_FitRec_vsPh_CB[i][j] = hfFitRecCB->makeTH2D(Form("FitRec %s vs Ph^{rec} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"#phi [deg]",Form("Fit - Rec %s",vartitle[j].c_str()),PhBins,VarBins.at(j),Form("hFitRec%svsPhCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_FitRec_vsPh_TAPS[i][j] = hfFitRecTAPS->makeTH2D(Form("FitRec %s vs Ph^{rec} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"#phi [deg]",Form("Fit - Rec %s",vartitle[j].c_str()),PhBins,VarBins.at(j),Form("hFitRec%svsPhTAPS%s",varname[j].c_str(),partname[i].c_str()),true);
            //-- TrueFit
            h_TrueFit_vsEk_CB[i][j] = hfTrueFitCB->makeTH2D(Form("TrueFit %s vs Ek^{true} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"Ek [MeV]",Form("True - Fit %s",vartitle[j].c_str()),EkBins,VarBins.at(j),Form("hTrueFit%svsEkCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueFit_vsEk_TAPS[i][j] = hfTrueFitTAPS->makeTH2D(Form("TrueFit %s vs Ek^{true} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"Ek [MeV]",Form("True - Fit %s",vartitle[j].c_str()),EkBins,VarBins.at(j),Form("hTrueFit%svsEkTAPS%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueFit_vsTh_CB[i][j] = hfTrueFitCB->makeTH2D(Form("TrueFit %s vs Th^{true} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"#theta [deg]",Form("True - Fit %s",vartitle[j].c_str()),ThCBBins,VarBins.at(j),Form("hTrueFit%svsThCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueFit_vsTh_TAPS[i][j] = hfTrueFitTAPS->makeTH2D(Form("TrueFit %s vs Th^{true} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"#theta [deg]",Form("True - Fit %s",vartitle[j].c_str()),ThTABins,VarBins.at(j),Form("hTrueFit%svsThTAPS%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueFit_vsPh_CB[i][j] = hfTrueFitCB->makeTH2D(Form("TrueFit %s vs Ph^{true} CB %s",vartitle[j].c_str(),parttitle[i].c_str()),"#phi [deg]",Form("True - Fit %s",vartitle[j].c_str()),PhBins,VarBins.at(j),Form("hTrueFit%svsPhCB%s",varname[j].c_str(),partname[i].c_str()),true);
            h_TrueFit_vsPh_TAPS[i][j] = hfTrueFitTAPS->makeTH2D(Form("TrueFit %s vs Ph^{true} TAPS %s",vartitle[j].c_str(),parttitle[i].c_str()),"#phi [deg]",Form("True - Fit %s",vartitle[j].c_str()),PhBins,VarBins.at(j),Form("hTrueFit%svsPhTAPS%s",varname[j].c_str(),partname[i].c_str()),true);

        }
    }
    //--- The Ebeam and z-vertex histos
    h_TrueFit_zvert = hfTrueFit->makeTH2D("True-Fit z-vertex vs True","z_{true}","z_{true} - z_{fit}",BinSettings(50,-15,15),BinSettings(51,-10,10),"h_TrueFit_zvert",true);
    h_TrueRec_Ebeam = hfTrueRec->makeTH2D("TrueRec E_{beam} vs E_{beam}^{rec}","E_{beam}^{rec}","True - Rec E_{beam}",BinSettings(1600,0.,1600.),BinSettings(101,-10,10),"h_TrueRec_Ebeam",true);
    h_FitRec_Ebeam = hfFitRec->makeTH2D("FitRec E_{beam} vs E_{beam}^{rec}","E_{beam}^{rec}","Fit - Rec E_{beam}",BinSettings(1600,0.,1600.),BinSettings(101,-10,10),"h_FitRec_Ebeam",true);
    h_TrueFit_Ebeam = hfTrueFit->makeTH2D("TrueFit E_{beam} vs E_{beam}^{true}","E_{beam}^{true}","True - Fit E_{beam}",BinSettings(1600,0.,1600.),BinSettings(101,-10,10),"h_TrueFit_Ebeam",true);
    //--- The IM and MM histos
    h_IMeeg_True = hfOverview->makeTH1D("IM(eeg) True","IM(e^{+}e^{-}#gamma) [MeV]","",BinSettings(500,0.,1000.),"h_IMeeg_True",true);
    h_IMeeg_Rec = hfOverview->makeTH1D("IM(eeg) Rec","IM(e^{+}e^{-}#gamma) [MeV]","",BinSettings(500,0.,1000.),"h_IMeeg_Rec",true);
    h_IMeeg_Fit = hfOverview->makeTH1D("IM(eeg) Fit","IM(e^{+}e^{-}#gamma) [MeV]","",BinSettings(500,0.,1000.),"h_IMeeg_Fit",true);
    h_IMgg_True = hfOverview->makeTH1D("IM(gg) True","IM(#gamma#gamma) [MeV]","",BinSettings(500,0.,1000.),"h_IMgg_True",true);
    h_IMgg_Rec = hfOverview->makeTH1D("IM(gg) Rec","IM(#gamma#gamma) [MeV]","",BinSettings(500,0.,1000.),"h_IMgg_Rec",true);
    h_IMgg_Fit = hfOverview->makeTH1D("IM(gg) Fit","IM(#gamma#gamma) [MeV]","",BinSettings(500,0.,1000.),"h_IMgg_Fit",true);

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
        truerecpair.push_back(truepart);
        if(match){
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

void scratch_lheijken_checkkinfit::DoFitComparisons(const TParticlePtr fitprot, const TParticlePtr recprot, const TParticleList fitphot, const TParticleList recphot, const std::vector<std::vector<TParticlePtr>> matchtruerecopart, const std::vector<TParticlePtr> truepart, double zvertfit, double tw)
{
    double kinvar_r[3], kinvar_t[3], kinvar_fr[3], kinvar_tf[3];
    for(int i=0; i<3; i++){kinvar_r[i]=-1000; kinvar_t[i]=-1000; kinvar_fr[i]=-1000; kinvar_tf[i]=-1000;}
    vec3 vertshift{0,0,zvertfit};

    //-- compare proton, if there is a reconstructed proton
    if(recprot){
        //-- calculate fit - rec (for theta: shift the fit with the extracted vertex)
        kinvar_r[0] = recprot->Ek(); kinvar_r[1] = recprot->Theta()*radtodeg; kinvar_r[2] = recprot->Phi()*radtodeg;
        kinvar_fr[0] = fitprot->Ek() - kinvar_r[0];
        kinvar_fr[1] = ((fitprot->p + vertshift).Theta()*radtodeg - kinvar_r[1]);
        kinvar_fr[2] = fitprot->Phi()*radtodeg - kinvar_r[2];
        //-- fill the fit-rec histograms
        if(recprot->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            for(int i=0; i<nrKinVars;i++){
                h_FitRec_vsEk_CB[en_p][i]->Fill(kinvar_r[0],kinvar_fr[i],tw);
                h_FitRec_vsTh_CB[en_p][i]->Fill(kinvar_r[1],kinvar_fr[i],tw);
                h_FitRec_vsPh_CB[en_p][i]->Fill(kinvar_r[2],kinvar_fr[i],tw);
            }
        }
        else{
            for(int i=0; i<nrKinVars;i++){
                h_FitRec_vsEk_TAPS[en_p][i]->Fill(kinvar_r[0],kinvar_fr[i],tw);
                h_FitRec_vsTh_TAPS[en_p][i]->Fill(kinvar_r[1],kinvar_fr[i],tw);
                h_FitRec_vsPh_TAPS[en_p][i]->Fill(kinvar_r[2],kinvar_fr[i],tw);
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
        kinvar_tf[0] = kinvar_t[0] - fitprot->Ek();
        kinvar_tf[1] = kinvar_t[1] - fitprot->Theta()*radtodeg;
        kinvar_tf[2] = kinvar_t[2] - fitprot->Phi()*radtodeg;
        //--- fill the true-fit histograms
        if(fitprot->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
            for(int i=0; i<nrKinVars;i++){
                h_TrueFit_vsEk_CB[en_p][i]->Fill(kinvar_t[0],kinvar_tf[i],tw);
                h_TrueFit_vsTh_CB[en_p][i]->Fill(kinvar_t[1],kinvar_tf[i],tw);
                h_TrueFit_vsPh_CB[en_p][i]->Fill(kinvar_t[2],kinvar_tf[i],tw);
            }
        }
        else{
            for(int i=0; i<nrKinVars;i++){
                h_TrueFit_vsEk_TAPS[en_p][i]->Fill(kinvar_t[0],kinvar_tf[i],tw);
                h_TrueFit_vsTh_TAPS[en_p][i]->Fill(kinvar_t[1],kinvar_tf[i],tw);
                h_TrueFit_vsPh_TAPS[en_p][i]->Fill(kinvar_t[2],kinvar_tf[i],tw);
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
                kinvar_fr[0] = fitphot.at(i)->Ek() - kinvar_r[0];
                kinvar_fr[1] = ((fitphot.at(i)->p + vertshift).Theta()*radtodeg - kinvar_r[1]);
                kinvar_fr[2] = fitphot.at(i)->Phi()*radtodeg - kinvar_r[2];
                //-- calculate the true - fit
                kinvar_t[0] = matchtruepart->Ek(); kinvar_t[1] = matchtruepart->Theta()*radtodeg; kinvar_t[2] = matchtruepart->Phi()*radtodeg;
                kinvar_tf[0] = kinvar_t[0] - fitphot.at(i)->Ek();
                kinvar_tf[1] = kinvar_t[1] - fitphot.at(i)->Theta()*radtodeg;
                kinvar_tf[2] = kinvar_t[2] - fitphot.at(i)->Phi()*radtodeg;
                //-- fill the fit-rec histograms
                if(recphot.at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                    for(int j=0; j<nrKinVars;j++){
                        h_FitRec_vsEk_CB[ptype][j]->Fill(kinvar_r[0],kinvar_fr[j],tw);
                        h_FitRec_vsTh_CB[ptype][j]->Fill(kinvar_r[1],kinvar_fr[j],tw);
                        h_FitRec_vsPh_CB[ptype][j]->Fill(kinvar_r[2],kinvar_fr[j],tw);
                    }
                }
                else{
                    for(int j=0; j<nrKinVars;j++){
                        h_FitRec_vsEk_TAPS[ptype][j]->Fill(kinvar_r[0],kinvar_fr[j],tw);
                        h_FitRec_vsTh_TAPS[ptype][j]->Fill(kinvar_r[1],kinvar_fr[j],tw);
                        h_FitRec_vsPh_TAPS[ptype][j]->Fill(kinvar_r[2],kinvar_fr[j],tw);
                    }
                }
                //-- fill true-fit histograms
                if(fitphot.at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                    for(int j=0; j<nrKinVars;j++){
                        h_TrueFit_vsEk_CB[ptype][j]->Fill(kinvar_t[0],kinvar_tf[j],tw);
                        h_TrueFit_vsTh_CB[ptype][j]->Fill(kinvar_t[1],kinvar_tf[j],tw);
                        h_TrueFit_vsPh_CB[ptype][j]->Fill(kinvar_t[2],kinvar_tf[j],tw);
                    }
                }
                else{
                    for(int j=0; j<nrKinVars;j++){
                        h_TrueFit_vsEk_TAPS[ptype][j]->Fill(kinvar_t[0],kinvar_tf[j],tw);
                        h_TrueFit_vsTh_TAPS[ptype][j]->Fill(kinvar_t[1],kinvar_tf[j],tw);
                        h_TrueFit_vsPh_TAPS[ptype][j]->Fill(kinvar_t[2],kinvar_tf[j],tw);
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
            kinvar_fr[0] = fitphot.at(i)->Ek() - kinvar_r[0];
            kinvar_fr[1] = ((fitphot.at(i)->p + vertshift).Theta()*radtodeg - kinvar_r[1]);
            kinvar_fr[2] = fitphot.at(i)->Phi()*radtodeg - kinvar_r[2];
            //-- fill the fit-rec histograms
            if(recphot.at(i)->Candidate->FindCaloCluster()->DetectorType == Detector_t::Type_t::CB){
                for(int j=0; j<nrKinVars;j++){
                    h_FitRec_vsEk_CB[en_g][j]->Fill(kinvar_r[0],kinvar_fr[j],tw);
                    h_FitRec_vsTh_CB[en_g][j]->Fill(kinvar_r[1],kinvar_fr[j],tw);
                    h_FitRec_vsPh_CB[en_g][j]->Fill(kinvar_r[2],kinvar_fr[j],tw);
                }
            }
            else{
                for(int j=0; j<nrKinVars;j++){
                    h_FitRec_vsEk_TAPS[en_g][j]->Fill(kinvar_r[0],kinvar_fr[j],tw);
                    h_FitRec_vsTh_TAPS[en_g][j]->Fill(kinvar_r[1],kinvar_fr[j],tw);
                    h_FitRec_vsPh_TAPS[en_g][j]->Fill(kinvar_r[2],kinvar_fr[j],tw);
                }
            }
        }
    }
}

AUTO_REGISTER_PHYSICS(scratch_lheijken_checkkinfit)
