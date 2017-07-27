#pragma once

#include "analysis/input/DataReader.h"

#include "base/WrapTTree.h"
#include "base/types.h"

#include <string>
#include <set>

namespace ant {

class WrapTFileInput;
struct TEventData;

namespace expconfig {
namespace detector {
struct Trigger;
}}

namespace analysis {
namespace input {

class GoatReader : public DataReader {
protected:

    std::shared_ptr<const WrapTFileInput> inputfiles;
    std::shared_ptr<const expconfig::detector::Trigger> trigger;

    using trees_t = std::set<std::reference_wrapper<TTree>>;

    struct treeDetectorHitInput_t {
        struct treeHits_t : WrapTTree {
            treeHits_t(const std::string& branchNamePrefix) :
                WrapTTree(branchNamePrefix) {}
            ADD_BRANCH_T(ROOTArray<Int_t>,    Hits)
        };

        struct treeEnergyTime_t : treeHits_t {
            treeEnergyTime_t(const std::string& branchNamePrefix) :
                treeHits_t(branchNamePrefix) {}
            ADD_BRANCH_T(ROOTArray<Double_t>, Energy)
            ADD_BRANCH_T(ROOTArray<Double_t>, Time)
        };

        struct treeCluster_t : treeEnergyTime_t {
            treeCluster_t(const std::string& branchNamePrefix) :
                treeEnergyTime_t(branchNamePrefix) {}
            ADD_BRANCH_T(ROOTArray<Int_t>,    Cluster)
        };

        treeHits_t       MWPC{"MWPC"};
        treeEnergyTime_t PID{"PID"};
        treeEnergyTime_t Veto{"Veto"};
        treeCluster_t    NaI{"NaI"};
        treeCluster_t    BaF2{"BaF2"};

        bool LinkBranches(const WrapTFileInput& input, trees_t& trees);
        void Copy(TEventData& recon);
    };

    struct treeTaggerInput_t {
        struct tree_t : WrapTTree {
            ADD_BRANCH_T(ROOTArray<Int_t>,    taggedChannel)
            ADD_BRANCH_T(ROOTArray<Double_t>, taggedTime)
            ADD_BRANCH_T(ROOTArray<Double_t>, taggedEnergy)
        };
        tree_t t;

        bool LinkBranches(const WrapTFileInput& input, trees_t& trees);
        void Copy(TEventData& recon);
    };

    struct treeTriggerInput_t {
        struct tree_t : WrapTTree {
            ADD_BRANCH_T(Double_t, 	energySum)
            ADD_BRANCH_T(Int_t,     multiplicity)
            ADD_BRANCH_T(ROOTArray<Int_t>, triggerPattern)
            ADD_BRANCH_T(ROOTArray<Int_t>, errorModuleID)
            ADD_BRANCH_T(ROOTArray<Int_t>, errorModuleIndex)
            ADD_BRANCH_T(ROOTArray<Int_t>, errorCode)
            ADD_BRANCH_OPT_T(Bool_t,   helicity)
            ADD_BRANCH_OPT_T(Long64_t, MC_evt_id)
            ADD_BRANCH_OPT_T(Long64_t, MC_rnd_id)
        };

        struct treeEventParameters_t : WrapTTree {
            ADD_BRANCH_T(Int_t, eventNumber)
            ADD_BRANCH_T(Int_t, nReconstructed)
        };

        tree_t t;
        treeEventParameters_t tEventParams;

        bool LinkBranches(const WrapTFileInput& input, trees_t& trees);
        void Copy(TEventData& recon);
    };

    struct treeTrackInput_t {
        struct tree_t : WrapTTree {
            ADD_BRANCH_T(ROOTArray<Double_t>, clusterEnergy)
            ADD_BRANCH_T(ROOTArray<Double_t>, theta)
            ADD_BRANCH_T(ROOTArray<Double_t>, phi)
            ADD_BRANCH_T(ROOTArray<Double_t>, time)
            ADD_BRANCH_T(ROOTArray<Int_t>, clusterSize)
            ADD_BRANCH_T(ROOTArray<Int_t>, centralCrystal)
            ADD_BRANCH_T(ROOTArray<Int_t>, centralVeto)
            ADD_BRANCH_T(ROOTArray<Int_t>, detectors)
            //Charged detector energies
            ADD_BRANCH_T(ROOTArray<Double_t>, vetoEnergy)
            ADD_BRANCH_T(ROOTArray<Double_t>, MWPC0Energy)
            ADD_BRANCH_T(ROOTArray<Double_t>, MWPC1Energy)
            //TAPS PSA Short-gate Energy
            ADD_BRANCH_T(ROOTArray<Double_t>, shortEnergy)
            //Pseudo vertex information
            ADD_BRANCH_T(ROOTArray<Double_t>, pseudoVertexX)
            ADD_BRANCH_T(ROOTArray<Double_t>, pseudoVertexY)
            ADD_BRANCH_T(ROOTArray<Double_t>, pseudoVertexZ)
        };

        tree_t t;

        bool LinkBranches(const WrapTFileInput& input, trees_t& trees);
        void Copy(TEventData& recon);
    };

    treeDetectorHitInput_t treeDetectorHitInput;
    treeTaggerInput_t      treeTaggerInput;
    treeTriggerInput_t     treeTriggerInput;
    treeTrackInput_t       treeTrackInput;

    template<typename... Args>
    static void insert_trees(trees_t& trees, const Args&... args) {
        trees.insert({std::ref(static_cast<TTree&>(*args.Tree))...});
    }

    trees_t trees;
    long long current_entry;
    long long max_entries;
    bool init;

public:
    GoatReader(const std::shared_ptr<const WrapTFileInput>& rootfiles);
    virtual ~GoatReader();

    virtual reader_flags_t GetFlags() const override;
    virtual bool ReadNextEvent(event_t& event) override;

    double PercentDone() const override;
};

}
}
}
