#include <iostream>
#include <string>
#include <sstream>
#include <vector>


#include "simulation/mc/A2Cocktail.h"

#include "expconfig/ExpConfig.h"

#include "base/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "base/detail/tclap/ValuesConstraintExtra.h"



using namespace std;
using namespace ant;
using namespace ant::simulation::mc;

int main( int argc, char** argv )
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-cocktail - Pluto-based cocktail generator for A2-Physics", ' ', "0.1");

    auto cmd_outfile    = cmd.add<TCLAP::ValueArg<string>> ("o", "outfile",       "Data output file"            , true, "", "filename");
    auto cmd_numEvents  = cmd.add<TCLAP::ValueArg<int>>    ("n", "numEvents",     "Number of generated events"  , true, 0, "# events");

    auto cmd_numEnergyBins = cmd.add<TCLAP::ValueArg<int>>    ("N", "numEnergyBins", "Number of tagged photon energy bins", false, 0, "# energy bins");
    auto cmd_Emin          = cmd.add<TCLAP::ValueArg<double>> ("",  "Emin",          "Minimum taggeg photon energy"       , false, 0, "MeV");
    auto cmd_Emax          = cmd.add<TCLAP::ValueArg<double>> ("",  "Emax",          "Maximum taggeg photon energy"       , false, 0, "MeV");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Setup to determine tagged photon energy bins",false,"",&allowedsetupnames);

    auto cmd_noBulk     = cmd.add<TCLAP::SwitchArg>        ("",  "no-bulk",       "disable Pluto-Bulk-Interface", false);
    auto cmd_noUnstable = cmd.add<TCLAP::SwitchArg>        ("",  "no-unstable",   "don't save unstable particles", false);

    auto cmd_dataFiles  = cmd.add<TCLAP::MultiArg<string>> ("",  "datafiles",     "Xsection-data-files",          false,       "inputfiles");

    cmd.parse(argc, argv);

    string outfile_clean(cmd_outfile->getValue());
    if(ant::std_ext::string_ends_with(outfile_clean, ".root")) {
        outfile_clean = outfile_clean.substr(0,outfile_clean.size()-5);
    }

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


    A2Cocktail cocktail(outfile_clean,
                        energies,
                        !cmd_noUnstable->isSet(), !cmd_noBulk->isSet(),
                        cmd_dataFiles->getValue()
                        );

    auto nErrors = cocktail.Sample(cmd_numEvents->getValue());

    if(nErrors>0)
        LOG(WARNING) << "Events with error: " <<  nErrors;
    cocktail.Finish();
}


