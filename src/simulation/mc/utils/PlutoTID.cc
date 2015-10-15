#include "PlutoTID.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "tree/TDataRecord.h"

#include "TTree.h"
#include "TRandom2.h"

#include <ctime>

using namespace std;
using namespace ant;
using namespace ant::simulation::mc::utils;

const string PlutoTID::tidtree_name("data_tid");
const string PlutoTID::geant_tidtree_name("h12_tid");

template <typename T>
bool CheckExists(WrapTFile& file, const std::string& name) {
    T* obj = nullptr;
    return file.GetObject(name, obj);
}

void PlutoTID::AddTID(const std::string &filename)
{
    const auto random_bits = 4;

    WrapTFileOutput file(filename, WrapTFileOutput::mode_t::update, true);

    if(CheckExists<TTree>(file, tidtree_name)) {
        LOG(WARNING) << "TID tree already exists in " << filename;
        return;
    }

    TTree* data = nullptr;
    if(file.GetObject("data", data)) {


        TTree* data_tid = file.CreateInside<TTree>(tidtree_name.c_str(),"Ant-TID for pluto data");

        TID tid(std::time(nullptr), 0, {TID::Flags_t::MC});

        data_tid->Branch("tid",&tid);

        auto nEvents = data->GetEntries();

        if(nEvents >= 1 << (sizeof(TID::Lower)*8 - random_bits)) {
           throw std::runtime_error("Too many entries to fit into TID together with random bits.");
        }

        TRandom2 rng;
        rng.SetSeed();

        for(decltype(nEvents) i=0; i<nEvents; ++i) {

            unsigned r = floor(rng.Uniform(1 << random_bits));
            tid.Lower = (i << random_bits) + r;

            data_tid->Fill();
        }
    } else {
        LOG(WARNING) << "No pluto data Tree found";
    }

}

void PlutoTID::CopyTIDPlutoGeant(const string& pluto_filename, const string& geant_filename)
{
    WrapTFileInput input(pluto_filename);

    TTree* intree = nullptr;
    input.GetObject(tidtree_name,intree);

    if(intree) {
        WrapTFileOutput outfile(geant_filename,WrapTFileOutput::mode_t::update, true);
        intree->SetBranchStatus("tid",1);
        auto outtree = intree->CloneTree(-1,"fast SortBasketsByBranch");
        outtree->SetName(geant_tidtree_name.c_str());
        outtree->Write();

    } else {
        LOG(ERROR) << "TTree \"" << tidtree_name << "\" not found";
    }
}
