#pragma once

#include "analysis/physics/Physics.h"
#include <memory>

class TCanvas;

namespace ant {

class TH2CB;

namespace analysis {
namespace physics {

class XMasCB : public Physics {
protected:

    std::unique_ptr<TH2CB> hist = nullptr;
    std::unique_ptr<TH2CB> grid = nullptr;
    TCanvas* c = nullptr;
    unsigned n = 0;
    unsigned m = 0;

    unsigned skipevents = 0;
    std::string ext;
    const int w_px;


public:
    XMasCB(const std::string& name, PhysOptPtr opts);
    virtual ~XMasCB();

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
};

}
}
}
