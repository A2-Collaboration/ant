#pragma once

#include "physics/Physics.h"
#include "base/interval.h"
#include <map>
#include <string>
#include "TLorentzVector.h"
class TH1;
class TTree;
class TClonesArray;

namespace ant {
namespace analysis {
namespace physics {

class DeltaTreeGenerator: public Physics {
private:
    TTree* photonTree;
    double taggerEnergy;
    TClonesArray* reconstructed;
    TClonesArray* mctrue;
    TH1D*  taggerHits;
    TH1D*  mcgamma;
    TH1D*  recgamma;

public:
    DeltaTreeGenerator(const std::string& name, OptionsPtr opts);
    virtual ~DeltaTreeGenerator() {}

    // Physics interface
    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void Finish() override;
    virtual void ShowResult() override;
};

}
}
}
