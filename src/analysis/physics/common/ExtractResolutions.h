#pragma once

#include "analysis/physics/Physics.h"

#include "TVector3.h"

class TH2D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

//struct ValueSigma_t {
//    double value = 0.0;
//    double sigma = 0.0;
//    ValueSigma_t(const double v=0.0, const double s=0.0):
//        value(v), sigma(s) {}

//    ValueSigma_t(const ValueSigma_t&) = default;
//    ValueSigma_t(ValueSigma_t&&) = default;
//    ValueSigma_t& operator= (const ValueSigma_t&) = default;
//    ValueSigma_t& operator= (      ValueSigma_t&&) = default;

//};

//struct Resolutions {
//    ValueSigma_t E = {};
//    ValueSigma_t Theta = {};
//    ValueSigma_t Phi = {};

//    Resolutions() {}

//    Resolutions(const Resolutions&) = default;
//    Resolutions(Resolutions&&) = default;
//    Resolutions& operator= (const Resolutions&) = default;
//    Resolutions& operator= (      Resolutions&&) = default;

//};

class ExtractResolutions : public Physics {
public:

    double   b_DE      = {};
    double   b_DTheta  = {};
    double   b_DPhi    = {};
    unsigned b_Element = {};
    double   b_trueE   = {};
    TVector3 b_trueDir = {};
    TVector3 b_recDir  = {};

    TTree* tree = nullptr;

    Detector_t::Type_t det = Detector_t::Type_t::CB;

    ExtractResolutions(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}
}
}
