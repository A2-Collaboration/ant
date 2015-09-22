#include "PlutoTID.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "tree/TDataRecord.h"

#include "TTree.h"

using namespace std;
using namespace ant;
using namespace ant::simulation::mc::utils;

const string PlutoTID::tidtree_name("data_tid");

template <typename T>
bool CheckExists(WrapTFile& file, const std::string& name) {
    T* obj = nullptr;
    return file.GetObject(name, obj);
}

void PlutoTID::AddTID(const std::string &filename)
{
    WrapTFileOutput file(filename, WrapTFileOutput::mode_t::update, true);

    if(CheckExists<TTree>(file, tidtree_name)) {
        LOG(WARNING) << "TID tree already exists in " << filename;
        return;
    }

    TTree* data = nullptr;
    if(file.GetObject("data", data)) {


        TTree* data_tid = file.CreateInside<TTree>(tidtree_name.c_str(),"Ant-TID for pluto data");

        TID tid(0, 0, true);

        data_tid->Branch("tid",&tid);

        auto nEvents = data->GetEntries();

        for(decltype(nEvents) i=0; i<nEvents; ++i) {
            ++tid;
            data_tid->Fill();
        }
    } else {
        LOG(WARNING) << "No pluto data Tree found";
    }

}
