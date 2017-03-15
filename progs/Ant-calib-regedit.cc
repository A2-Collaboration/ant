#include "base/Logger.h"
#include "tclap/CmdLine.h"

#include "base/std_ext/system.h"
#include "base/std_ext/string.h"
#include "base/WrapTFile.h"
#include "tree/TCalibrationData.h"
#include "tree/TAntHeader.h"

#include "TTree.h"

#include <iostream>

using namespace ant;
using namespace ant::std_ext;
using namespace std;

void convert(string setupfolder);
void build_index(const std::vector<string>& files);

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("Ant-calib-regedit - manage calib database on-disk format", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");

    TCLAP::SwitchArg cmd_mode_convert("","convert","Convert old database in given setup folders (needs already created structure)", false);
    TCLAP::SwitchArg cmd_mode_index("","index","Build Data File TID Ranges index from Ant files", false);

    cmd.xorAdd(cmd_mode_convert, cmd_mode_index);

    auto cmd_givenstrings  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");

    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    if(cmd_mode_convert.isSet()) {
        if(cmd_givenstrings->getValue().empty())
            LOG(ERROR) << "No input folders given!";

        for(auto setupfolder : cmd_givenstrings->getValue())
            convert(setupfolder);
    }

    if(cmd_mode_index.isSet()) {
        build_index(cmd_givenstrings->getValue());
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

void build_index(const std::vector<string>& files) {
    WrapTFileOutput outf("index.root");

    TTree* tree = outf.CreateInside<TTree>("index", "Run file index");

    TID start;
    TID end;
    string b_file;

    tree->Branch("start", &start);
    tree->Branch("end",   &end);
    tree->Branch("file",  &b_file);


    for(const auto& file : files) {

        WrapTFileInput f(file);

        TAntHeader* header(nullptr);
        f.GetObject("AntHeader", header);

        if(header) {
            start = header->FirstID;
            end   = header->LastID;
            b_file = file;
            tree->Fill();
            LOG(INFO) << "Added " << file << ": " << start << " - " << end;
        } else {
            LOG(ERROR) << "No AntHeader in " << file;
        }
    }
}
