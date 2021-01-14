#pragma once

#include "physics/Physics.h"
#include "base/interval.h"
#include "base/piecewise_interval.h"
#include "utils/Matcher.h"

#include <string>

namespace ant {
namespace analysis {
namespace physics {

/**
 * @brief A class for plotting Ek True-Rec vs Rec for
 * different particle types (currently only p, e+, e- and g).
 * Can be used with single particle MC or normal pluto MC
 */

class TrueRecCheck_ClusterE : public Physics {

protected:
    const interval<double> CBThetaWindow;
    const interval<double> TAPSThetaWindow;
    const PiecewiseInterval<double> CBHemisphereGap;

    static const int nrParticles = 4;
    enum parttype{en_p=0, en_ep, en_em, en_g};
    enum dettype{en_cb=0, en_ta};
    static constexpr auto radtodeg = std_ext::radian_to_degree(1.0);

    TH2D *h_Ek_TrueRecvsRec[2][nrParticles];
    TH1D *h_PairedOpAngle[2][nrParticles];

public:
    TrueRecCheck_ClusterE(const std::string& name, OptionsPtr opts);
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;

};


}}}
