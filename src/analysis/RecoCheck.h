#ifndef RECOCHECK_H
#define RECOCHECK_H

#include "AntPhysics.h"
#include "base/interval.h"

class TH1D;
class TH2D;

namespace ant {
namespace analysis {

class RecoCheck: public Physics {
protected:
    TH1D* angle_diff;
    TH1D* n_unmatched;

    interval<radian_t> cb_angle;

public:
    RecoCheck();
    virtual ~RecoCheck() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif
