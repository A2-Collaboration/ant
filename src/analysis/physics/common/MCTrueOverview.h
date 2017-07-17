#pragma once

#include "analysis/physics/Physics.h"
#include "base/ParticleTypeTree.h"
#include "tree/TParticle.h"
#include "base/WrapTTree.h"

#include <map>

namespace ant {
namespace analysis {
namespace physics {


class MCTrueOverview : public Physics {
protected:

    struct CBTAPS_Multiplicity {
        TH1D* h_CBTAPS = nullptr;
        CBTAPS_Multiplicity(const HistogramFactory& histFac, const int nParticles);
        void Fill(const TParticleTree_t& particles);
    };

    struct perChannel_t {
        perChannel_t(const HistogramFactory& histFac,
                     const TParticleTree_t& ptree);

        void Fill(const TParticleTree_t& ptree, const TCandidateList& cands);
        void Show(canvas& c) const;

    protected:

        struct histnode_t {
            using typeptr_t = const ParticleTypeDatabase::Type*;
            histnode_t(std::unique_ptr<const HistogramFactory> histFacPtr,
                       const std::vector<typeptr_t>& leafTypes);
            void Fill(const TParticle& p);
            void Show(canvas& c) const;

        protected:
            struct perType_t {
                perType_t(const HistogramFactory& HistFac,
                          const ParticleTypeDatabase::Type& type);
                TH2D* h_EkTheta = nullptr;
                TH2D* h_EkTheta_EtaPrime = nullptr;
            };

            std::map<typeptr_t, perType_t> hists;

        };

        using histtree_t = Tree<histnode_t>::node_t;
        histtree_t histtree;
        void traverse_tree_and_fill(const histtree_t& histtree, const TParticleTree_t& ptree) const;

        TH1D* h_CBEsum_true = nullptr;
        TH1D* h_CBEsum_rec = nullptr;

        CBTAPS_Multiplicity mult;

        struct tree_t : WrapTTree {
            ADD_BRANCH_T(std::vector<TLorentzVector>,  LV)
            ADD_BRANCH_T(std::vector<double>,          M)
            ADD_BRANCH_T(std::vector<double>,          Ek)
            ADD_BRANCH_T(std::vector<double>,          Theta)
            ADD_BRANCH_T(std::vector<double>,          Phi)
        };
        tree_t tree;
    };

    std::map<ParticleTypeTreeDatabase::Channel, std::unique_ptr<perChannel_t>> channels;

public:
    MCTrueOverview(const std::string& name, OptionsPtr opts);

    virtual void ProcessEvent(const TEvent& event, manager_t& manager) override;
    virtual void ShowResult() override;
};

}
}
}
