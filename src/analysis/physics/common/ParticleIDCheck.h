#pragma once

#include "analysis/physics/Physics.h"
#include "base/interval.h"

#include <string>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class ParticleIDCheck : public Physics {
protected:

    struct branch_hists {
        branch_hists(SmartHistFactory& HistFac,const std::string& name);
        TH1D* hist;
        void Fill(const TEventData& data);
    };

    branch_hists mctrue;
    branch_hists rec;

    std::vector< std::tuple<interval<double>,interval<double>,TH2D*> > bananas;

public:
    ParticleIDCheck(const std::string& name,OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
