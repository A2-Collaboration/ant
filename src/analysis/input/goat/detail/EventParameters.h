#pragma once

#include "InputModule.h"
#include "Rtypes.h"

#include <vector>

namespace ant {
namespace analysis {
namespace input {

struct EventParameters : BaseInputModule {
    Int_t EventNumber = 0;
    Int_t nReconstructed = 0;

    EventParameters();
    virtual ~EventParameters();

    bool SetupBranches(TreeRequestManager&& input_files);
    void GetEntry();

};
}
}
}
