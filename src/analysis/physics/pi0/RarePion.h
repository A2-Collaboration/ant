#pragma once

#include "physics/Physics.h"

#include <map>

class TH1D;

namespace ant {
namespace analysis {
namespace physics {

class RarePion: public Physics {
protected:
    TH1D* ggim;
    TH1D* gggim;
    TH1D* ggggim;
    TH1D* nphotons;
    TH1D* nphotons_anglecut;
    TH1D* im_3g;
    TH1D* minAngle;
    TH1D* minAngle2g;
    TH1D* minAngle3g;
    TH1D* splitoffangle;

    TH2D* theta_vs_En;


    //std::map<const ant::ParticleTypeDatabase::Type*, TH1D*> EHists;

public:
    RarePion(const std::string& name, OptionsPtr opts);
    virtual ~RarePion() {}

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
