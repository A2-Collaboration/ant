#include <iostream>
#include <string>
#include <sstream>
#include <vector>


#include "simulation/mc/A2Cocktail.h"

#include "base/CmdLine.h"
#include "base/std_ext/string.h"



using namespace std;
using namespace ant::simulation::mc;

int main( int argc, char** argv )
{

    TCLAP::CmdLine cmd("Ant-cocktail - Pluto-based cocktail generator for A2-Physics", ' ', "0.1");

    auto cmd_outfile    = cmd.add<TCLAP::ValueArg<string>> ("o", "outfile",       "Data output file"            , true,    "", "filename");
    auto cmd_numEvents  = cmd.add<TCLAP::ValueArg<int>>    ("N", "numEvents",     "Number of generated events"  , true,     0, "# events");
    auto cmd_numBins    = cmd.add<TCLAP::ValueArg<int>>    ("n", "numEnergyBins", "Number of taggerbins"        , false,   47, "# energy bins");
    auto cmd_Emin       = cmd.add<TCLAP::ValueArg<double>> ("",  "Emin",          "Minimum Tagger energy"       , false, 1420, "MeV");
    auto cmd_Emax       = cmd.add<TCLAP::ValueArg<double>> ("",  "Emax",          "Maximum Tagger energy"       , false, 1580, "MeV");

    auto cmd_noBulk     = cmd.add<TCLAP::SwitchArg>        ("",  "no-bulk",       "disable Pluto-Bulk-Interface", false);
    auto cmd_noUnstable = cmd.add<TCLAP::SwitchArg>        ("",  "no-unstable",   "don't save unstable particles", false);

    auto cmd_dataFiles  = cmd.add<TCLAP::MultiArg<string>> ("",  "datafiles",     "Xsection-data-files",          false,       "inputfiles");

    cmd.parse(argc, argv);

    string outfile_clean(cmd_outfile->getValue());
    if(ant::std_ext::string_ends_with(outfile_clean, ".root")) {
        outfile_clean = outfile_clean.substr(0,outfile_clean.size()-5);
    }

    A2Cocktail cocktail(outfile_clean,
                        cmd_Emin->getValue() / 1000.0, cmd_Emax->getValue() / 1000.0,  // Cocktailclass takes GeV
                        cmd_numBins->getValue(),
                        !cmd_noUnstable->isSet(), !cmd_noBulk->isSet(),
                        cmd_dataFiles->getValue()
                        );

    cout << "Events with error: " << cocktail.Sample(cmd_numEvents->getValue()) << endl;
    cocktail.Finish();
}


