#include "analysis/physics/pi0/Pi0Dalitz.h"
#include "analysis/plot/HistogramFactory.h"
#include "analysis/physics/Physics.h"
#include "analysis/utils/Combinatorics.h"

#include "base/Logger.h"
#include "base/WrapTFile.h"
#include "base/std_ext/system.h"

#include "tclap/CmdLine.h"

#include "tree/TAntHeader.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TRint.h"

#include <iostream>
#include <memory>

/**
 *  A program for creating the pi0 TFF from the output of the Pi0Dalitz physics class
 **/

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;


// global useful stuff
static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

static const int nrProtCases = 3;
enum prcases{en_allpr=0, en_wpr, en_nopr};
static const int nrPcuts = 5;
static const double Pcuts[] = {0.,0.01,0.02,0.05,0.1};
static const int nrPIDCases = 3;
enum pidcases{en_allpid=0, en_samepid, en_diffpid};
static const int nrIMeebins = 20;
static const double IMeebinedges[] = {0.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,100.,110.,120.,500.};

// functions
void CreateHist(const Pi0Dalitz::TFFTree_t &tfftree, const HistogramFactory hf);
void TreeHistos(const HistogramFactory hf);
int FindIMeeBin(double imee);

// histograms
TH1D *h_Probgg[nrProtCases], *h_Probeeg[nrProtCases];
TH1D *h_IMgg[nrProtCases][nrPcuts];
TH1D *h_IMeeg[nrProtCases][nrPcuts][nrPIDCases], *h_IMeeg_BinCut[nrProtCases][nrPcuts][nrPIDCases], *h_IMeeg_IMeeBins[nrProtCases][nrPcuts][nrPIDCases][nrIMeebins];
TH1D *h_IMee;

void TreeHistos(const HistogramFactory hf)
{
    //-- Specific bin settings
    vector<double> vIMeebinedges;
    for(int i=0; i<nrIMeebins+1; i++)
        vIMeebinedges.push_back(IMeebinedges[i]);
    auto BinsIMee = VarBinSettings(vIMeebinedges);
    auto BinsIMgg = BinSettings(500,-0.5,499.5);
    auto BinsIMeeg = BinSettings(500,-0.5,499.5);

    //-- Names and titles prefixes
    string prcasename[] = {"allprot","withprot","noprot"};
    string prcasetitle[] = {"all proton cases","with proton","no proton"};
    string pidcasename[] = {"allpid","samepid","diffpid"};
    string pidcasetitle[] = {"all pid cases","same pid","different pid"};

    //-- The main folders
    auto hfIMgg = new HistogramFactory("IMgg",hf,"");
    auto hfIMeeg = new HistogramFactory("IMeeg",hf,"");

    //-- The histograms for IMgg
    for(int i=0; i<nrProtCases; i++){
        auto hfProtCaseTemp = new HistogramFactory(prcasename[i],*hfIMgg,"");
        h_Probgg[i] = hfProtCaseTemp->makeTH1D(Form("P(#chi^{2}), 2#gamma decay, %s",prcasetitle[i].c_str()),"P(#chi^{2}_{KF})","",BinSettings(1000,0.,1.),Form("hProbgg_%s",prcasename[i].c_str()),true);
        for(int j=0; j<nrPcuts; j++){
            auto hfPCutsTemp = new HistogramFactory(Form("PCut_%03d",(int)(Pcuts[j]*100)),*hfProtCaseTemp,"");
            h_IMgg[i][j] = hfPCutsTemp->makeTH1D(Form("IM(#gamma#gamma), %s, P(#chi^{2})>%.2f",prcasetitle[i].c_str(),Pcuts[j]),"IM(#gamma#gamma)","",BinsIMgg,Form("hIMgg_%s_pc%03d",prcasename[i].c_str(),(int)(Pcuts[j]*100)),true);
            delete hfPCutsTemp;
        }
        delete hfProtCaseTemp;
    }
    //-- The histograms for IMeeg
    for(int i=0; i<nrProtCases; i++){
        auto hfProtCaseTemp = new HistogramFactory(prcasename[i],*hfIMeeg,"");
        h_Probeeg[i] = hfProtCaseTemp->makeTH1D(Form("P(#chi^{2}), e^{+}e^{-}#gamma decay, %s",prcasetitle[i].c_str()),"P(#chi^{2}_{KF})","",BinSettings(1000,0.,1.),Form("hProbeeg_%s",prcasename[i].c_str()),true);
        for(int j=0; j<nrPcuts; j++){
            auto hfPCutsTemp = new HistogramFactory(Form("PCut_%03d",(int)(Pcuts[j]*100)),*hfProtCaseTemp,"");
            for(int k=0; k<nrPIDCases; k++){
                auto hfPIDCaseTemp = new HistogramFactory(Form("PIDcase_%s",pidcasename[k].c_str()),*hfPCutsTemp,"");
                h_IMeeg[i][j][k] = hfPIDCaseTemp->makeTH1D(Form("IM(e^{+}e^{-}#gamma), %s, P(#chi^{2})>%.2f, %s",prcasetitle[i].c_str(),Pcuts[j],pidcasetitle[k].c_str()),"IM(e^{+}e^{-}#gamma)","",BinsIMeeg,Form("hIMeeg_%s_pc%03d_%s",prcasename[i].c_str(),(int)(Pcuts[j]*100),pidcasename[k].c_str()),true);
                h_IMeeg_BinCut[i][j][k] = hfPIDCaseTemp->makeTH1D(Form("IM(e^{+}e^{-}#gamma) %s, P(#chi^{2})>%.2f, %s, IM(ee)#in[%.1f,%.1f)",prcasetitle[i].c_str(),Pcuts[j],pidcasetitle[k].c_str(),vIMeebinedges.at(0),vIMeebinedges.at(nrIMeebins-1)),"IM(e^{+}e^{-}#gamma)","",BinsIMeeg,Form("hIMeeg_BinCut_%s_pc%03d_%s",prcasename[i].c_str(),(int)(Pcuts[j]*100),pidcasename[k].c_str()),true);
                for(int l=0; l<nrIMeebins; l++){
                    h_IMeeg_IMeeBins[i][j][k][l] = hfPIDCaseTemp->makeTH1D(Form("IM(e^{+}e^{-}#gamma), %s, P(#chi^{2})>%.2f, %s, IM(ee)#in[%.1f,%.1f)",prcasetitle[i].c_str(),Pcuts[j],pidcasetitle[k].c_str(),vIMeebinedges.at(l),vIMeebinedges.at(l+1)),"IM(e^{+}e^{-}#gamma)","",BinsIMeeg,Form("hIMeeg_%s_pc%03d_%s_bin%d",prcasename[i].c_str(),(int)(Pcuts[j]*100),pidcasename[k].c_str(),l),true);
                }
                delete hfPIDCaseTemp;
            }
            delete hfPCutsTemp;
        }
        delete hfProtCaseTemp;
    }
    //-- This is just a dummy IMee histogram to be used in later analysis since it stores the IMee bins used
    h_IMee = hfIMeeg->makeTH1D(Form("Dummy IM(e^{+}e^{-}), all proton cases, P(#chi^{2})>%.2f, all pid cases",Pcuts[1]),"IM(e^{+}e^{-})","",BinsIMee,"hIMee_dummy",true);

}

void CreateHist(const Pi0Dalitz::TFFTree_t& tfftree, const HistogramFactory hf)
{
    //-- Create the histograms
    TreeHistos(hf);

    //-- Loop over the tree entries (N.B. here there's one entry per event)
    const long long int max_entries = tfftree.Tree->GetEntries();
    long long entry = 0;
    while (entry < max_entries) {
        tfftree.Tree->GetEntry(entry++);
        double percdone = 100*(double)entry/(double)max_entries;
        //if(entry%1000000==0) std::cout<<percdone<<"% of events done."<<std::endl;

        //-- Fetch some of the tree info (easier access)
        bool IsDalDec = tfftree.TBIsDalDec;
        bool Is2gDec = tfftree.TBIs2gDec;
        bool HasProton = tfftree.TBHasProtonCand;
        vector<double> IMeegs = tfftree.TBIMeegs;
        vector<double> IMggs = tfftree.TBIMggs;
        vector<double> IMees = tfftree.TBIMees;
        vector<bool> samePIDs = tfftree.TBsamePIDs;
        vector<double> KFprobs = tfftree.TBKFprobs;
        vector<double> TWeights = tfftree.TBTaggWeights;

        int prcaseind = en_allpr;
        if(HasProton) prcaseind = en_wpr;
        else prcaseind = en_nopr;
        vector<int> prcaseinds {en_allpr,prcaseind};

        //-- Loop over the tagger hits
        for(int i=0; i<(int)TWeights.size(); i++){

            //-- gg events
            if(Is2gDec){
                for(int prc : prcaseinds){
                    h_Probgg[prc]->Fill(KFprobs.at(i),TWeights.at(i));

                    for(int j=0; j<nrPcuts; j++){
                        if(KFprobs.at(i) >= Pcuts[j]){
                            h_IMgg[prc][j]->Fill(IMggs.at(i),TWeights.at(i));
                        }
                    }
                }
            }

            //-- daldec events
            if(IsDalDec){
                int pidcaseind = en_allpid;
                if(samePIDs.at(i)) pidcaseind = en_samepid;
                else pidcaseind = en_diffpid;
                vector<int> pidcaseinds {en_allpid,pidcaseind};

                for(int prc : prcaseinds){
                    h_Probeeg[prc]->Fill(KFprobs.at(i),TWeights.at(i));
                    for(int j=0; j<nrPcuts; j++){
                        if(KFprobs.at(i) >= Pcuts[j]){
                            for(int pidc : pidcaseinds){
                                h_IMeeg[prc][j][pidc]->Fill(IMeegs.at(i),TWeights.at(i));
                                if(IMeebinedges[0] <= IMees.at(i) && IMeebinedges[nrIMeebins-1] > IMees.at(i))
                                    h_IMeeg_BinCut[prc][j][pidc]->Fill(IMeegs.at(i),TWeights.at(i));
                                int imeebin = FindIMeeBin(IMees.at(i));
                                if(imeebin>-1){
                                    h_IMeeg_IMeeBins[prc][j][pidc][imeebin]->Fill(IMeegs.at(i),TWeights.at(i));
                                }
                                //else {
                                //    if(prc==0 && j==0 && pidc==0)
                                //        LOG(INFO)<<"IM(ee) is outside range: "<<IMees.at(i)<<", PC="<<KFprobs.at(i);
                                //}
                            }
                        }
                    }
                }
                //-- fill the dummy IMee histogram
                if(KFprobs.at(i)>=Pcuts[1])
                    h_IMee->Fill(IMees.at(i),TWeights.at(i));
            }
        }
    }
}

int FindIMeeBin(double imee)
{
    for(int i=0; i<nrIMeebins; i++){
        if(IMeebinedges[i]<=imee && IMeebinedges[i+1]>imee)
            return i;
    }
    return -1;
}

int main( int argc, char** argv )
{
    SetupLogger();

    //-- Read from the command line
    TCLAP::CmdLine cmd("Create the pi0 TFF from the output of the Pi0Dalitz physics class", ' ', "0.1");
    auto cmd_batchmode = cmd.add<TCLAP::MultiSwitchArg>("b", "batch","Run in batch mode (no ROOT shell afterwards)", false);
    auto cmd_input = cmd.add<TCLAP::ValueArg<string>>("i","input","Input file",true,"","intputfiles");
    auto cmd_output = cmd.add<TCLAP::ValueArg<string>>("o", "output", "Output file", false, "", "filename");
    cmd.parse(argc, argv);

    //-- Create the output file
    argc=0; // prevent TRint to parse any cmdline
    TRint app("Pi0TFF",&argc,argv,nullptr,0,true);
    unique_ptr<WrapTFileOutput> masterFile;
    if(cmd_output->isSet()) {
        // cd into masterFile upon creation
        masterFile = std_ext::make_unique<WrapTFileOutput>(cmd_output->getValue(), true);
    }
    HistogramFactory hf("Pi0TFF");

    //-- Open the input file
    WrapTFileInput input(cmd_input->getValue());

    //-- Get the tree
    Pi0Dalitz::TFFTree_t tfftree;
    if (!input.GetObject("Pi0Dalitz/tff_variables", tfftree.Tree))
        runtime_error("Can't find the tree Pi0Dalitz/tff_variables!");
    tfftree.LinkBranches();
    //-- Creating histograms from the tree
    CreateHist(tfftree, hf);

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

