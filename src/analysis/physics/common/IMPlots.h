#pragma once

#include "analysis/physics/Physics.h"

namespace ant {
namespace analysis {
namespace physics {

class IMPlots : public Physics {
protected:
    struct hist_set {
        std::vector<TH1D*> m;
        void Fill(unsigned ngamma, double mm);
        hist_set(const std::string& pref, SmartHistFactory& hf, std::size_t n=8);
    };

    hist_set cb;
    hist_set taps;

public:
    IMPlots(PhysOptPtr opts);

    void ProcessEvent(const data::Event &event) override;
    void Finish() override;
    void ShowResult() override;
};

}
}
}
