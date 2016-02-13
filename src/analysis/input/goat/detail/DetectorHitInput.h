#pragma once

#include "InputModule.h"
#include "Rtypes.h"
#include <vector>
#include "TMath.h"

#define GTreeHits_MAX 1024

namespace ant {
namespace analysis {
namespace input {

struct DetectorHitInput : BaseInputModule {

    struct Item_t {
        Item_t(const std::string& name);
        const std::string Name;
        Int_t nHits;
        Int_t Hits[GTreeHits_MAX];
        Int_t Cluster[GTreeHits_MAX];
        Double_t Energy[GTreeHits_MAX];
        Double_t Time[GTreeHits_MAX];
        void SetupBranches(TTree* tree);
    };

    Item_t NaI;
    Item_t PID;
    Item_t MWPC;
    Item_t BaF2;
    Item_t Veto;

    DetectorHitInput();
    virtual ~DetectorHitInput();

    bool SetupBranches(TreeRequestManager&& input_files);
    void GetEntry();
};
}
}
}
