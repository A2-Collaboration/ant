#include <iostream>
#include <string>
#include <sstream>
#include <vector>


#include "simulation/mc/A2Channels.h"

#include "expconfig/ExpConfig.h"

#include "TH1D.h"
#include "TRint.h"
#include "TStyle.h"


#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"
#include "base/iterators.h"

#include "analysis/plot/root_draw.h"


using namespace std;
using namespace ant;
using namespace ant::simulation::mc;

int main( int argc, char** argv )
{
    SetupLogger();

    TCLAP::CmdLine cmd("show_cocktail_database - plot current cocktail database", ' ', "0.1");


    auto cmd_numEnergyBins = cmd.add<TCLAP::ValueArg<int>>    ("N", "numEnergyBins", "Number of tagged photon energy bins", false, 0, "# energy bins");
    auto cmd_Emin          = cmd.add<TCLAP::ValueArg<double>> ("",  "Emin",          "Minimum taggeg photon energy"       , false, 0, "MeV");
    auto cmd_Emax          = cmd.add<TCLAP::ValueArg<double>> ("",  "Emax",          "Maximum taggeg photon energy"       , false, 0, "MeV");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Setup to determine tagged photon energy bins",false,"",&allowedsetupnames);

    auto cmd_dataFiles  = cmd.add<TCLAP::MultiArg<string>> ("",  "datafiles",     "Xsection-data-files",          false,       "inputfiles");

    cmd.parse(argc, argv);

    vector<double> energies; // in GeV
    constexpr double MeVtoGeV = 1.0/1000.0;

    if(cmd_setup->isSet()) {

        if(cmd_numEnergyBins->isSet() || cmd_Emin->isSet() || cmd_Emax->isSet()) {
            LOG(ERROR) << "Either specify setup name or manually choose energy bins";
            return 1;
        }

        ExpConfig::Setup::ManualName = cmd_setup->getValue();
        auto setup = ExpConfig::Setup::GetLastFound();
        if(auto tagger = setup->GetDetector<TaggerDetector_t>()) {
            for(unsigned ch=0;ch<tagger->GetNChannels();ch++) {
                energies.push_back(tagger->GetPhotonEnergy(ch)*MeVtoGeV); // convert to GeV!
            }
        }
        else {
            LOG(ERROR) << "Specified setup did not provide tagger detector";
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
            energies.push_back(energy*MeVtoGeV); // convert to GeV!
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

    A2ChannelManager aman(cmd_dataFiles->getValue());
    auto channels = aman.GetChannels();
    vector<TH1D*> histograms;

    //production channels
    for (const auto& channel: channels)
    {
        TH1D* histptr = new TH1D(channel.c_str(),channel.c_str(),energies.size()-1,&energies[0]);
        for(int i = 1 ; i <= histptr->GetNbinsX() ; ++i)
            histptr->SetBinContent(i,aman.Xsection(channel,histptr->GetBinCenter(i)));
        histograms.push_back(histptr);
    }
    //total
    TH1D* totalptr = new TH1D("total","Total Cross-section",energies.size()-1,&energies[0]);
    for(int i = 1 ; i <= totalptr->GetNbinsX() ; ++i)
        totalptr->SetBinContent(i,aman.TotalXsection(totalptr->GetBinCenter(i)));
    histograms.push_back(totalptr);

    hstack sumplot("sum", "");
    auto cit = getCirculatIterator(ColorPalette::Colors.begin(), ColorPalette::Colors.end());
    for (auto histptr: histograms)
    {
        histptr->SetXTitle("E_{#gamma} [GeV]");
        histptr->SetYTitle("#sigma_{tot} [#mub]");
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


    canvas("total") << drawoption("nostack") << padoption::set(padoption_t::Legend) << padoption::set(padoption_t::LogY)
                    << sumplot << endc;

    app->Run(kTRUE);

}

