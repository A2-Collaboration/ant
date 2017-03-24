#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "mc/database/Query.h"

#include "expconfig/ExpConfig.h"

#include "TH1D.h"
#include "TRint.h"
#include "TStyle.h"


#include "tclap/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "tclap/ValuesConstraintExtra.h"
#include "base/std_ext/iterators.h"
#include "analysis/utils/ParticleTools.h"

#include "analysis/plot/RootDraw.h"

#include "root-addons/analysis_codes/hstack.h"

#include "detail/McAction.h"


using namespace std;
using namespace ant;
using namespace ant::mc::data;
using namespace ant::analysis;

int main( int argc, char** argv )
{
    SetupLogger();

    TCLAP::CmdLine cmd("show_cocktail_database - plot current cocktail database", ' ', "0.1");

    const auto allowedTargetNames = getAllowedTargetNames();
    TCLAP::ValuesConstraintExtra<std::vector<string>> allowedTargetsConstrain(allowedTargetNames);
    auto cmd_target = cmd.add<TCLAP::ValueArg<string>>("","target","Target Particle",false,"proton",&allowedTargetsConstrain);


    auto cmd_numEnergyBins = cmd.add<TCLAP::ValueArg<int>>    ("N", "numEnergyBins", "Number of tagged photon energy bins", false, 0, "# energy bins");
    auto cmd_Emin          = cmd.add<TCLAP::ValueArg<double>> ("",  "Emin",          "Minimum taggeg photon energy"       , false, 0, "MeV");
    auto cmd_Emax          = cmd.add<TCLAP::ValueArg<double>> ("",  "Emax",          "Maximum taggeg photon energy"       , false, 0, "MeV");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Setup to determine tagged photon energy bins",false,"",&allowedsetupnames);

    auto cmd_dataFiles  = cmd.add<TCLAP::MultiArg<string>> ("",  "datafiles",     "Xsection-data-files",          false,       "inputfiles");

    cmd.parse(argc, argv);



    vector<double> energies;


    if(cmd_setup->isSet()) {

        if(cmd_numEnergyBins->isSet() || cmd_Emin->isSet() || cmd_Emax->isSet()) {
            LOG(ERROR) << "Either specify setup name or manually choose energy bins";
            return 1;
        }

        const auto& setupname = cmd_setup->getValue();

        try {
            ExpConfig::Setup::SetByName(setupname);
            auto tagger = ExpConfig::Setup::GetDetector<TaggerDetector_t>();
            for(unsigned ch=0;ch<tagger->GetNChannels();ch++) {
                energies.push_back(tagger->GetPhotonEnergy(ch));
            }
        }
        catch(ExpConfig::Exception e) {
            LOG(ERROR) << "Specified setup '" << setupname << "' did not provide a tagger: " << e.what();
            return 1;
        }
    }
    else {
        const double Emin = cmd_Emin->getValue();
        const double Emax = cmd_Emax->getValue();
        const double numBins = cmd_numEnergyBins->getValue();
        const double dE = (Emax-Emin)/numBins;
        for(unsigned bin=0;bin<numBins;bin++) {
            const double energy = Emin + (bin+0.5)*dE;
            energies.push_back(energy);
        }
    }

    if(energies.empty()) {
        LOG(ERROR) << "Please specify some photon energy bins, either manually or by setup name";
        return 1;
    }

    sort(energies.begin(),energies.end());

    int fake_argc=1;
    char* fake_argv[2];
    fake_argv[0] = argv[0];
    auto app = new TRint("Ant-calib",&fake_argc,fake_argv);

    auto selector = mc::data::Query::GetSelector(allowedTargets.at(cmd_target->getValue()));


    auto channels = Query::GetProductionChannels(selector);
    vector<TH1D*> histograms;

    //production channels
    for (const auto& channel: channels)
    {
        auto chTree = ParticleTypeTreeDatabase::Get(channel);
        auto chString = utils::ParticleTools::GetDecayString(chTree);
        auto chShortString = utils::ParticleTools::GetDecayString(chTree,false);

        TH1D* histptr = new TH1D(chShortString.c_str(),chString.c_str(),energies.size()-1,&energies[0]);
        for(int i = 1 ; i <= histptr->GetNbinsX() ; ++i)
        {
            auto bc = histptr->GetBinCenter(i);
            auto q = Query::Xsection(channel,bc);
            histptr->SetBinContent(i,q);
        }
        histograms.push_back(histptr);
    }
    //total
    TH1D* totalptr = new TH1D("total","Total Cross-section",energies.size()-1,&energies[0]);
    for(int i = 1 ; i <= totalptr->GetNbinsX() ; ++i)
    {
        auto bc = totalptr->GetBinCenter(i);
        auto q = Query::TotalXsection(bc,selector);
        totalptr->SetBinContent(i,q);
    }
    histograms.push_back(totalptr);

    hstack sumplot("sum", "");
    auto cit = std_ext::getCircularIterator(ColorPalette::Colors.begin(), ColorPalette::Colors.end());
    for (auto histptr: histograms)
    {
        histptr->GetXaxis()->SetNdivisions(4);
        histptr->GetYaxis()->SetNdivisions(3);
        histptr->SetXTitle("E_{#gamma} [MeV]");
        histptr->SetYTitle("#sigma_{tot} [#mub]");
        histptr->SetLabelSize(0.08,"X");
        histptr->SetLabelSize(0.08,"Y");
        histptr->SetMarkerStyle(7);
        histptr->SetLineColor(*cit);
        sumplot << histptr;
        cit.next();
    }

    gStyle->SetOptStat(kFALSE);
    canvas c("overview");
    c << drawoption("LP");// << drawoption("same");
    //c << totalptr;
    for (auto hist: histograms)
        c << hist;
    c << endc;

    //WTF u not WORK!!!!!????   ROOT stinkt!
//    sumplot.GetXaxis()->SetTitle("E_{#gamma} [MeV]");
//    sumplot.GetYaxis()->SetTitle("#sigma_{tot} [#mub]");

    canvas("total") << drawoption("nostack") // << padoption::Legend << padoption::LogY
                    << &sumplot << endc;

    app->Run(kTRUE);

}

