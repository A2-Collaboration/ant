/**
  * @file Ant-addTID.cc
  * @brief Program to add TID tree to a pluto ROOT file or a geant ROOT file.
  *
  * Adds a TTree to the given pluto ROOT file, containing a branch of TID with one ID for each pluto event (entry in the "data" TTree).
  *
  * The "copy-to-geant" mode copies a TID TTree from a pluto ROOT file to a geant ROOT file.
  * Useful if the TID tree was not added after the pluto file was created and geant has aleady simulated the pluto data.
  * This should be used with caution and if you are very sure that the pluto and geant file belong together.
  */

#include "base/Logger.h"
#include "tclap/CmdLine.h"
#include "base/WrapTFile.h"
#include "base/std_ext/memory.h"

#include "TH1.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

using namespace std;
using namespace ant;

enum class HistMismatch_t {
    Match,
    Entries,
    Axes,
    BinContent
};

bool compare(const TH1* a, const TH1* b) {

    bool result = true;
    stringstream message;

    if(a->GetEntries() != b->GetEntries()) {
        result = false;
        message << " Entries ";
    }

    if(a->GetNbinsX() != b->GetNbinsX()) {
        result = false;
        message << " NXbins";
    }

    if(a->GetXaxis()->GetXmin() != b->GetXaxis()->GetXmin()) {
        result = false;
        message << " Xmin";
    }

    if(a->GetXaxis()->GetXmax() != b->GetXaxis()->GetXmax()) {
        result = false;
        message << " Xmax";
    }

    unsigned nbins = 0;
    for(int i=0; i<=a->GetNbinsX()+1; ++i) {
        if(a->GetBinContent(i) != b->GetBinContent(i)) {
            nbins++;
        }
    }

    if(nbins != 0) {
        result = false;
        message << " Bin Content(" << nbins << ")";
    }

    if(result == false)
        LOG(ERROR) << "Mismatch: " << a->GetName() << ":"  << message.str();

    return result;
}

int main (int argc, char** argv)
{
    SetupLogger();

    TCLAP::CmdLine cmd("Ant-compareHists - Compare histograms bin by bin", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");

    auto cmd_h1 = cmd.add<TCLAP::MultiArg<string>>("","h1","Name of 1D Histogram to compare",true,"");

    // unlabeled multi arg must be the last element added, and interprets everything as a input file
    auto cmd_inputfiles  = cmd.add<TCLAP::UnlabeledMultiArg<string>>("inputfiles","Ant files with histograms",true,"inputfiles");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    vector<unique_ptr<WrapTFileInput>> files;

    for(const auto& filename : cmd_inputfiles->getValue()) {
        try {
            auto file = std_ext::make_unique<WrapTFileInput>(filename);
            files.emplace_back(move(file));
        } catch (std::exception& e) {
            LOG(ERROR) << e.what();
            return EXIT_FAILURE;
        }
    }

    for(const auto& hname : cmd_h1->getValue()) {
        for(auto i=files.begin(); i!=files.end();++i) {
            for(auto j=next(i); j!=files.end();++j) {

                TH1* h1=nullptr;
                (*i)->GetObject(hname, h1);

                if(!h1) {
                    LOG(ERROR) << hname << " not found in " << cmd_inputfiles->getValue().at(size_t(distance(files.begin(),i)));
                    return EXIT_FAILURE;
                }

                TH1* h2=nullptr;
                (*j)->GetObject(hname, h2);
                if(!h2) {
                    LOG(ERROR) << hname << " not found in " << cmd_inputfiles->getValue().at(size_t(distance(files.begin(),j)));
                    return EXIT_FAILURE;
                }

                if(!compare(h1,h2)) {
                    LOG(ERROR) << hname << " mismatch between "
                               << cmd_inputfiles->getValue().at(size_t(distance(files.begin(),i)))
                               << " and "
                               << cmd_inputfiles->getValue().at(size_t(distance(files.begin(),j)));
                    return EXIT_FAILURE;
                }
            }
        }
    }

    return EXIT_SUCCESS;

}
