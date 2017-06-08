#pragma once
#include "analysis/physics/Physics.h"
#include "base/Tree.h"
#include "base/interval.h"
#include "base/std_ext/math.h"
#include "base/WrapTTree.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "Rtypes.h"
#include "tree/TSimpleParticle.h"
#include <map>


class TH1D;
class TH2D;
class TH3D;

namespace ant {

namespace analysis {
namespace physics {



class Pi0EtaLostG : public Physics {
public:
    struct Tree_t : WrapTTree {
        Tree_t();

        ADD_BRANCH_T(std::vector<TLorentzVector>, lostV, 5)
        ADD_BRANCH_T(std::vector<int>, lostid, 5)
    };

    struct lostCounter {
        TH2D* lostPhotons;
        TH2D* producedPhotons;
        TH2D* normalized = nullptr;
        const std::string prefix;

        lostCounter(HistogramFactory& hf, const std::string& name);
        lostCounter(lostCounter&&) = default;
        lostCounter& operator=(lostCounter&&) = default;

        void Finish(HistogramFactory& hf);
        void FillLost(const TParticle& p);
        void FillProduced(const TParticle& p);

    };

    lostCounter lostPhtons;
    lostCounter lostProtons;
    std::map<const ParticleTypeDatabase::Type*, lostCounter> counters;


protected:
    TTree*  tree;
    Tree_t t;

public:

    Pi0EtaLostG(const std::string& name, OptionsPtr opts);
    virtual ~Pi0EtaLostG();

    void ProcessEvent(const TEvent& event, manager_t&) override;
    void ShowResult() override;
    void Finish() override;
};

}
}
}
