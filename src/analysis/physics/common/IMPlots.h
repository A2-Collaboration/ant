#pragma once

#include "analysis/physics/Physics.h"
class TH1D;
class TTree;

namespace ant {
namespace analysis {
namespace physics {

class IMPlots : public Physics {
public:
    struct hist_set {
        std::vector<TH1D*> m;
        void Fill(unsigned ngamma, double mm);
        hist_set(const std::string& pref, SmartHistFactory& hf, std::size_t n=8);
        unsigned MinNGamma() const noexcept { return 2;}
        unsigned MaxNGamma() const noexcept { return m.size()+2; }
    };

//    hist_set cb;
//    hist_set taps;
    hist_set all;

public:
    IMPlots(const std::string& name, PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};


class Symmetric2Gamma : public Physics {
protected:
    TH1D* h_symIM = nullptr;
    TTree* tree   = nullptr;

    double b_IM;
    double b_E;
    double b_E1;
    double b_E2;
    double b_theta1;
    double b_theta2;
    double b_phi1;
    double b_phi2;

    double perc = 0.2;

public:
    Symmetric2Gamma(const std::string& name, PhysOptPtr opts);
    virtual ~Symmetric2Gamma();

    void ProcessEvent(const data::Event &event) override;
    void ShowResult() override;
};

}
}
}
