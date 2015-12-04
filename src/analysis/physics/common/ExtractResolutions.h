#pragma once

#include "analysis/physics/Physics.h"

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

    TH2D* E_offset      = nullptr;
    TH2D* E_sigma       = nullptr;
    TH2D* Theta_offset  = nullptr;
    TH2D* Theta_sigma   = nullptr;
    TH2D* Phi_offset    = nullptr;
    TH2D* Phi_sigma     = nullptr;
    TTree* tree         = nullptr;

    double b_E = 0;
    double b_Theta = 0;
    double b_Phi = 0;
    unsigned b_Element = 0;

    Detector_t::Type_t det = Detector_t::Type_t::CB;

    ExtractResolutions(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

}
}
}
