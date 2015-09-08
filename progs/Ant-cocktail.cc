#include <iostream>
#include <string>
#include <sstream>
#include <vector>


#include "simulation/mc/A2Cocktail.h"

#include "base/CmdLine.h"



using namespace std;
using namespace ant::simulation::mc;

int main( int argc, char** argv )
{

    TCLAP::CmdLine cmd("Ant-cocktail - Pluto-based cocktail generator for A2-Physics", ' ', "0.1");

    auto cmd_outfile    = cmd.add<TCLAP::ValueArg<string>> ("o", "outfile",       "Data output file"            , true,    "", "filename");
    auto cmd_numEvents  = cmd.add<TCLAP::ValueArg<int>>    ("N", "numEvents",     "Number of generated events"  , true,     0, "# events");
    auto cmd_numBins    = cmd.add<TCLAP::ValueArg<int>>    ("n", "numEnergyBins", "Number of taggerbins"        , false,   48, "# energy bins");
    auto cmd_Emin       = cmd.add<TCLAP::ValueArg<double>> ("",  "Emin",          "Minimum Tagger energy"       , false, 1.42, "GeV");
    auto cmd_Emax       = cmd.add<TCLAP::ValueArg<double>> ("",  "Emax",          "Maximum Tagger energy"       , false, 1.58, "GeV");

    auto cmd_noBulk     = cmd.add<TCLAP::SwitchArg>        ("",  "no-bulk",       "disable Pluto-Bulk-Interface", false);
    auto cmd_noUnstable = cmd.add<TCLAP::SwitchArg>        ("",  "no-unstable",   "disable Pluto-Bulk-Interface", false);

    auto cmd_dataFiles  = cmd.add<TCLAP::MultiArg<string>> ("",  "datafiles",     "Xsection-data-files",          false,       "inputfiles");

    cmd.parse(argc, argv);

    A2Cocktail cocktail(cmd_outfile->getValue(),
                        cmd_Emin->getValue(), cmd_Emax->getValue(),
                        cmd_numBins->getValue(),
                        !cmd_noUnstable->isSet(), !cmd_noBulk->isSet(),
                        cmd_dataFiles->getValue()
                        );

    cout << "Events with error: " << cocktail.Sample(cmd_numEvents->getValue()) << endl;
    cocktail.Finish();
}


