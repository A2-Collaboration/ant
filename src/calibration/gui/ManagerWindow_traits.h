#pragma once

#include "Manager_traits.h"

namespace ant {
namespace calibration {
namespace gui {

class ManagerWindowGUI_traits : public ManagerWindow_traits {
public:
    struct Mode {
        Mode() :
            gotoNextSlice(true),
            autoContinue(true),
            autoFinish(false),
            showEachFit(true),
            skipStoreFit(false),
            channelStep(1),
            requestChannel(-1)
        {}

        bool gotoNextSlice;
        bool autoContinue;
        bool autoFinish;
        bool showEachFit;
        bool skipStoreFit;
        int  channelStep;
        int  requestChannel;
    };

    virtual Mode& GetMode() =0;
    virtual void SetProgressMax(unsigned slices, unsigned channels) =0;
    virtual void SetProgress(unsigned slice, unsigned channel) =0;
    virtual void SetFinishMode(bool flag) =0;

protected:
    ~ManagerWindowGUI_traits() = default;
};


}
}
}
