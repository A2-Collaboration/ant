#include <iostream>
#include <string>
#include <sstream>

#include "tclap/CmdLine.h"
#include "base/std_ext/string.h"
#include "base/Logger.h"
#include "tclap/ValuesConstraintExtra.h"

#include "TFile.h"

using namespace std;

int main( int argc, char** argv )
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-repair-keys - try to repair files with broken keys", ' ', "0.1");

    auto cmd_infile    = cmd.add<TCLAP::ValueArg<string>> ("i", "infile",       "root file input"          , true, "", "filename");

    cmd.parse(argc, argv);

    string filename(cmd_infile->getValue());

    TFile* f = new TFile(filename.c_str(),"UPDATE");

    f->Write();
    f->Close();

    return EXIT_SUCCESS;
}
