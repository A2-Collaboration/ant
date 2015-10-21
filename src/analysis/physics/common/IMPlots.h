#pragma once

#include "analysis/physics/Physics.h"

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

}
}
}
