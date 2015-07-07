#ifndef SPLITCHECK_H
#define SPLITCHECK_H

#include "AntPhysics.h"
#include "base/interval.h"
#include <map>
#include <vector>

class TH1D;
class TH2D;

namespace ant {
namespace analysis {

class SplitCheck: public Physics {
protected:
    std::map<int, std::vector<TH1D*>> gamma_rank;
    TH1D* big_small_angle;
    TH2D* small_theta_phi;
    TH2D* big_theta_phi;

    TH2D* IMsmall_theta_phi;
    TH2D* IMbig_theta_phi;

public:
    SplitCheck();
    virtual ~SplitCheck() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif
