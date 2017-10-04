#include <iostream>
#include <string>
#include <sstream>
#include <vector>


#include "mc/pluto/PlutoGenerator.h"
#include "mc/pluto/utils/PlutoTID.h"

#include "expconfig/ExpConfig.h"

#include "tclap/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "tclap/ValuesConstraintExtra.h"

#include "detail/McAction.h"

using namespace std;
using namespace ant;
using namespace ant::mc::pluto;

int main( int argc, char** argv )
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-cocktail - Pluto-based cocktail generator for A2-Physics", ' ', "0.1");

    auto cmd_outfile    = cmd.add<TCLAP::ValueArg<string>> ("o", "outfile",       "Data output file"          , true, "", "filename");
    auto cmd_numEvents  = cmd.add<TCLAP::ValueArg<int>>    ("n", "numEvents",     "Number of generated events", true, 0,  "# events");

    auto cmd_numEnergyBins = cmd.add<TCLAP::ValueArg<int>>    ("N", "numEnergyBins", "Number of tagged photon energy bins", false, 0, "# energy bins");
    auto cmd_Emin          = cmd.add<TCLAP::ValueArg<double>> ("",  "Emin",          "Minimum tagged photon energy",        false, 0.1, "MeV");
    auto cmd_Emax          = cmd.add<TCLAP::ValueArg<double>> ("",  "Emax",          "Maximum tagged photon energy",        false, 0.1, "MeV");

    TCLAP::ValuesConstraintExtra<decltype(ExpConfig::Setup::GetNames())> allowedsetupnames(ExpConfig::Setup::GetNames());
    auto cmd_setup = cmd.add<TCLAP::ValueArg<string>>("s","setup","Setup to determine tagged photon energy bins",false,"",&allowedsetupnames);

    const auto allowedTargetNames = getAllowedTargetNames();
    TCLAP::ValuesConstraintExtra<std::vector<string>> allowedTargetsConstrain(allowedTargetNames);
    auto cmd_target = cmd.add<TCLAP::ValueArg<string>>("","target","Target Particle",false,"proton",&allowedTargetsConstrain);

    auto cmd_noBulk     = cmd.add<TCLAP::SwitchArg>        ("",  "no-bulk",       "disable Pluto-Bulk-Interface",  false);
    auto cmd_noUnstable = cmd.add<TCLAP::SwitchArg>        ("",  "no-unstable",   "don't save unstable particles", false);

    auto cmd_noTID      = cmd.add<TCLAP::SwitchArg>        ("",  "noTID",   "Don't add TID tree for the events",   false);
    auto cmd_verbose    = cmd.add<TCLAP::ValueArg<int>>    ("v", "verbose", "Verbosity level (0..9)",              false, 0, "int");

    cmd.parse(argc, argv);

    string outfile(cmd_outfile->getValue());

    vector<double> energies; // in MeV

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
        if(!cmd_numEnergyBins->isSet()) {
            LOG(ERROR) << "Please specify some photon energy bins, either manually or by setup name";
            return 1;
        }
        const double Emin = cmd_Emin->getValue();
        const double Emax = cmd_Emax->getValue();
        if(Emin > Emax) {
            LOG(ERROR) << "Specified maximum energy is greater than the minimum energy";
            return 1;
        }
        const double numBins = cmd_numEnergyBins->getValue();
        const double dE = (Emax-Emin)/numBins;
        for(unsigned bin=0;bin<numBins;bin++) {
            const double energy = Emin + (bin+0.5)*dE;
            energies.push_back(energy); // convert to GeV!
        }
    }

    if(energies.empty()) {
        LOG(ERROR) << "Energy bins empty. Please specify a useful energy range and bins either manually or by setup name";
        return 1;
    }

    // scope that the Cocktail output file is properly closed before adding TID tree
    {
        auto selector = mc::data::Query::GetSelector(allowedTargets.at(cmd_target->getValue()));
        Cocktail cocktail(outfile,
                          energies,
                          !cmd_noUnstable->isSet(),
                          !cmd_noBulk->isSet(),
                          cmd_verbose->getValue(),
                          "1.0 / x",
                          selector);

        auto nErrors = cocktail.Sample(cmd_numEvents->getValue());

        if(nErrors>0)
            LOG(WARNING) << "Events with error: " <<  nErrors;
    }

    // add TID tree for the generated events
    if(!cmd_noTID->isSet()) {
        LOG(INFO) << "Add TID tree to the output file";
        mc::pluto::utils::PlutoTID::AddTID(outfile);
    }

    return EXIT_SUCCESS;
}
