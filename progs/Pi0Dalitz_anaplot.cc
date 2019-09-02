#include "analysis/physics/pi0/Pi0Dalitz.h"

#include "analysis/plot/HistogramFactory.h"
#include "analysis/physics/Physics.h"

#include "base/Logger.h"
#include "base/WrapTFile.h"
#include "base/std_ext/system.h"

#include "tclap/CmdLine.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TRint.h"

#include <iostream>
#include <memory>

/**
 * Analyse and plot the output of the
 * analysis class of the pi0->e+e-g channel
 **/

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

// functions
void CreateHistos(const HistogramFactory hf);
void DoTrueMCStuff(const int WhichMC, const TLorentzVector &ep, const TLorentzVector &em, const vector<TLorentzVector> &g);
void DoRecoCandStuff(const int cut, const TParticleList &recocand, const double& tw);

// histograms
TH1D *h_IMeegTrue, *h_IsPi0eegMC;
const int nrCuts = 2;
TH1D *h_IMeegReco[nrCuts], *h_IMggReco[nrCuts];

static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);


int main( int argc, char** argv )
{
    SetupLogger();

    //-- Read from the command line
    TCLAP::CmdLine cmd("Analysing and plotting the output of Pi0Dalitz analysis class", ' ', "0.1");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b", "batch","Run in batch mode (no ROOT shell afterwards)", false);
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","intputfiles");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o", "output", "Output file", false, "", "filename");
    cmd.parse(argc, argv);

    //-- Create the output file and the histograms to be stored there
    argc=0; // prevent TRint to parse any cmdline
    TRint app("Pi0Dalitz_anaplot",&argc,argv,nullptr,0,true);
    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }
    HistogramFactory hf("Pi0Dalitz_anaplot");
    CreateHistos(hf);

    //-- Open the input file and get the tree
    WrapTFileInput input(cmd_input->getValue());
    Pi0Dalitz::AnalysTree_t anatree;
    if (!input.GetObject("Pi0Dalitz/analysis_variables", anatree.Tree))
        runtime_error("Can't find the tree Pi0Dalitz/analysis_variables!");
    anatree.LinkBranches();

    //-- Loop over the tree entries (N.B. here there's one entry per event)
    const long long int max_entries = anatree.Tree->GetEntries();
    long long entry = 0;
    while (entry < max_entries) {
        anatree.Tree->GetEntry(entry++);

        //--- Fetch some of the tree info (easier access)
        vector<double> Tweights = anatree.TBPrRndWeight;
        vector<double> CaCaloE = anatree.TBCanCaloE;
        vector<double> CaThe = anatree.TBCanTheta;
        vector<double> CaPhi = anatree.TBCanPhi;
        vector<bool> CaInCB = anatree.TBCanInCB;
        vector<TLorentzVector> TrueVecGammas = anatree.TBTrueVecGammas;

        //--- Do some stuff with the true MC info
        if(anatree.TBMCpi0eeg){
            DoTrueMCStuff(1, anatree.TBTrueVecep, anatree.TBTrueVecem, TrueVecGammas);
        }

        //--- Create LorentzVectors for each candidate assuming no mass
        TParticleList RecoParticles;
        int nrCanInCB = 0;
        for(unsigned i=0; i<CaCaloE.size(); i++){
            TParticle p(ParticleTypeDatabase::Photon,CaCaloE.at(i),CaThe.at(i),CaPhi.at(i));
            RecoParticles.emplace_back(make_shared<TParticle>(ParticleTypeDatabase::Photon,CaCaloE.at(i),CaThe.at(i),CaPhi.at(i)));
            if(CaInCB.at(i)){
                nrCanInCB++;
            }
        }

        //--- Loop over the tagger hits
        for(unsigned i=0; i<Tweights.size(); i++){

            if(CaCaloE.size()==3 && nrCanInCB==3){
                DoRecoCandStuff(0,RecoParticles,Tweights.at(i));
            }

            if(CaCaloE.size()==4 && (nrCanInCB==3 || nrCanInCB==4))
                DoRecoCandStuff(1,RecoParticles,Tweights.at(i));
        }
    }

    if(!cmd_batchmode->isSet()){
        if(!std_ext::system::isInteractive()){
            LOG(INFO) << "No TTY attached. Not starting ROOT shell.";
        }
        else{
            if(masterFile)
                LOG(INFO) << "Close ROOT properly to write data to disk.";
            app.Run(kTRUE); // really important to return...
            if(masterFile)
                LOG(INFO) << "Writing output file...";
            masterFile = nullptr;   // and to destroy the master WrapTFile before TRint is destroyed
        }
    }

    return 0;
}

void CreateHistos(const HistogramFactory hf)
{
    auto hfTrueMC = new HistogramFactory("TrueMC",hf,"");
    auto hfCandChecks = new HistogramFactory("CandChecks",hf,"");

    h_IsPi0eegMC = hfTrueMC->makeTH1D("Is Pi0eeg","","",BinSettings(10,-0.5,9.5),"h_IsPi0eeg",true);
    h_IMeegTrue = hfTrueMC->makeTH1D("IMeeg True","IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),"h_IMeegTrue",true);

    string cutname[] = {"3Cand","4Cand"};
    string cuttitle[] = {"3 candidates","4 candidates"};
    for(int i=0; i<nrCuts; i++){
        h_IMeegReco[i] = hfCandChecks->makeTH1D(Form("IMeeg %s",cuttitle[i].c_str()),"IM(e^{+}e^{-}#gamma)","",BinSettings(250,0.,1000.),Form("h_IMeeg%s",cutname[i].c_str()),true);
        h_IMggReco[i] = hfCandChecks->makeTH1D(Form("IMgg %s",cuttitle[i].c_str()),"IM(#gamma#gamma)","",BinSettings(250,0.,1000.),Form("h_IMgg%s",cutname[i].c_str()),true);
    }
}

void DoTrueMCStuff(const int WhichMC, const TLorentzVector &ep, const TLorentzVector &em, const vector<TLorentzVector> &g)
{
    h_IsPi0eegMC->Fill(WhichMC);
    if(g.size()>0)
        h_IMeegTrue->Fill((ep+em+g.at(0)).M());
}

void DoRecoCandStuff(const int cut, const TParticleList &recocand, const double &tw)
{
    if(cut==0) {
        h_IMeegReco[cut]->Fill((*recocand.at(0)+*recocand.at(1)+*recocand.at(2)).M(),tw);
        h_IMggReco[cut]->Fill((*recocand.at(0)+*recocand.at(1)).M(),tw);
        h_IMggReco[cut]->Fill((*recocand.at(1)+*recocand.at(2)).M(),tw);
        h_IMggReco[cut]->Fill((*recocand.at(2)+*recocand.at(0)).M(),tw);
    }
    if(cut==1){
        h_IMeegReco[cut]->Fill((*recocand.at(0)+*recocand.at(1)+*recocand.at(2)).M(),tw);
        h_IMeegReco[cut]->Fill((*recocand.at(1)+*recocand.at(2)+*recocand.at(3)).M(),tw);
        h_IMeegReco[cut]->Fill((*recocand.at(2)+*recocand.at(3)+*recocand.at(0)).M(),tw);
    }
}

