#include "base/Logger.h"
#include "base/CmdLine.h"

#include "base/std_ext/system.h"
#include "base/std_ext/string.h"

#include <iostream>

using namespace ant;
using namespace std;

void create(vector<string> inputfiles);
void convert(string setupfolder);

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-calib-regedit - manage calib database on-disk format", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");

    TCLAP::SwitchArg cmd_mode_create("","create","Create the database skeleton from given files", false);
    TCLAP::SwitchArg cmd_mode_convert("","convert","Convert old database in given setup folders (needs already created structure)", false);
    cmd.xorAdd(cmd_mode_convert, cmd_mode_create);

    auto cmd_givenstrings  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(cmd_mode_create.isSet())
        create(cmd_givenstrings->getValue());
    else if(cmd_mode_convert.isSet()) {
        for(auto setupfolder : cmd_givenstrings->getValue())
            convert(setupfolder);
    }
}

void convert(string setupfolder) {
    LOG(INFO) << "Converting " << setupfolder;
    for(auto treefile : std_ext::system::lsFiles(setupfolder+"/calibration")) {
        if(!std_ext::string_ends_with(treefile, ".root"))
            continue;
        LOG(INFO) << "Opening " << treefile;
    }
}

void create(vector<string> inputfiles) {
    // ensure that nothing is there already
}
