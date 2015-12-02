#include "base/Logger.h"
#include "base/CmdLine.h"

#include "base/std_ext/system.h"
#include "base/std_ext/string.h"
#include "base/WrapTFile.h"
#include "tree/TCalibrationData.h"

#include "TTree.h"

#include <iostream>

using namespace ant;
using namespace ant::std_ext;
using namespace std;

void convert(string setupfolder);

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-calib-regedit - manage calib database on-disk format", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");

    auto cmd_mode_convert = cmd.add<TCLAP::SwitchArg>("","convert","Convert old database in given setup folders (needs already created structure)", false);

    auto cmd_givenstrings  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(cmd_mode_convert->isSet()) {
        if(cmd_givenstrings->getValue().empty())
            LOG(ERROR) << "No input folders given!";

        for(auto setupfolder : cmd_givenstrings->getValue())
            convert(setupfolder);
    }
}

string OutFilename(const Long64_t n) {
    return formatter() << setw(4) << setfill('0') << n << ".root";
}

void SplitTree(TTree* tree, const string& path) {


    TCalibrationData* data = nullptr;

    tree->SetBranchAddress("cdata", &data);

    Long64_t n = 0;

    for( n=0; n<tree->GetEntries(); ++n) {

        const string current_file(path+"/"+OutFilename(n));
        WrapTFileOutput of(current_file);

        tree->GetEntry(n);

        auto nbytes = of.WriteObject(data, "cdata");
        LOG(INFO) << "Written " << nbytes << " calibration data to file";
    }

    --n;

    if(n >= 0) {
        system::exec(formatter() << "cd " << path << " && ln -s -T -f " << OutFilename(n) << " current");
    }
}

void convert(string setupfolder) {
    LOG(INFO) << "Converting " << setupfolder;
    const string base_dir(setupfolder+"/calibration");

    for(auto treefile : std_ext::system::lsFiles(base_dir, ".root")) {

        LOG(INFO) << "Opening " << treefile;

        WrapTFileInput infile(treefile);

        auto trees = infile.GetListOf<TTree>();

        for(auto tree : trees) {
            string name = tree->GetName();
            LOG(INFO) << "Found TTree " << name;
            if(string_starts_with(name, "calibration-")) {
                auto calibname = name;
                removesubstr(calibname, "calibration-");
                LOG(INFO) << "Calibration: " << calibname;

                std_ext::replace(calibname, "/", "_");

                const string path = base_dir + "/" + calibname + "/DataDefault";

                system::exec(formatter() << "mkdir -p " << path );
                LOG(INFO) << "mkdir " << path;

                SplitTree(tree, path);

            }
        }

    }
}
