#ifndef BASIC_H
#define BASIC_H

#include "AntPhysics.h"
#include <vector>
#include <map>

class TH1D;
class TH2D;

namespace ant {
namespace analysis {

class Basic: public Physics {
protected:
    TH2D* banana;
    TH1D* particles;
    TH1D* tagger;
    TH1D* ntagged;

    TH1D* cbesum;

    std::map<int, std::map<int, TH1D*>> nGammaImEvent;


    std::map<const ParticleTypeDatabase::Type*, TH1D*> numParticleType;

public:
    Basic(const mev_t energy_scale=1000.0);
    virtual ~Basic() {}
    void ProcessEvent(const Event &event);
    void Finish();
    void ShowResult();
};

}
}
#endif
