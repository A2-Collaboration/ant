#pragma once

#include "TTree.h"

namespace ant {

/**
 * @brief The WrapTTree struct simplifies TTree handling in physics classes
 *
 * Example usage (see also etaprime_omega_gamma.h):
 *
 * struct treeTest_t : WrapTTree {
 *   ADD_BRANCH_T(bool,   IsSignal)
 *   ADD_BRANCH_T(double, KinFitChi2)
 *   // ... more branches using macro ADD_BRANCH_T
 * };
 *
 * treeTest_t treeTest;
 *
 * ... in constructor, setup for "writing" ...
 *
 *    treeTest.AddBranches(HistFac.makeTTree("test"));
 *
 * ... get/set branches via operator() ...
 *
 *    treeTest.IsSignal() = b_IsSignal;
 *
 * ... fill (or get entry) ...
 *
 *   treeTest.Tree->Fill();
 *
 */
struct WrapTTree {
    /**
     * @brief Tree to be used as usual TTree
     */
    TTree* Tree = nullptr;

    /**
     * @brief CreateBranches prepares the instance for filling the TTree
     * @param tree the tree to be filled
     */
    void CreateBranches(TTree* tree) {
        Tree = tree;
        // little trick to access the protected method
        struct TTree_trick : TTree {
            using TTree::BranchImpRef;
        };
        auto tree_trick = (TTree_trick*)Tree;
        for(ROOT_branch_t b : branches) {
            // dereference ValuePtr here to pointer to value
            tree_trick->BranchImpRef(b.Name.c_str(), b.ROOTClass, b.ROOTType, *b.ValuePtr, 32000, 99);
        }
    }

    /**
     * @brief LinkBranches prepares the instance for reading the TTree
     * @param tree the tree to read from
     */
    void LinkBranches(TTree* tree) {
        Tree = tree;
        for(ROOT_branch_t b : branches) {
            // copied from TTree::SetBranchAddress<T>
            // but handling of classes should be tested better...
            if(b.ROOTClass)
                Tree->SetBranchAddress(b.Name.c_str(),b.ValuePtr,0,b.ROOTClass,b.ROOTType,true);
            else
                Tree->SetBranchAddress(b.Name.c_str(),*b.ValuePtr,0,b.ROOTClass,b.ROOTType,false);
        }
    }

    template<typename T>
    struct Branch_t {
        Branch_t(WrapTTree* wraptree, const std::string& name) :
            Name(name),
            Value(new T())
        { wraptree->branches.emplace_back(ROOT_branch_t::Make(Name, std::addressof(Value))); }
        ~Branch_t() { delete Value; }
        Branch_t(const Branch_t&) = delete;
        Branch_t(Branch_t&&) = delete;
        Branch_t& operator= (const Branch_t&) = delete;
        Branch_t& operator= (Branch_t&&) = delete;

        const std::string Name;
        T* Value;

        operator T& () { return *Value; }
        operator const T& () const { return *Value; }
        T& operator= (const T& v) { *Value = v; return *Value; }
        T& operator= (T&& v) { *Value = v; return *Value; }
        // if you need to call methods of T, sometimes operator() is handy
        T& operator() () { return *Value; }
        const T& operator() () const { return *Value; }
    };

protected:

    template<typename T>
    friend struct Branch_t;

    struct ROOT_branch_t {
        template<typename T>
        static ROOT_branch_t Make(const std::string& name, T** valuePtr) {
            ROOT_branch_t b;
            b.Name = name;
            b.ValuePtr = (void**)valuePtr;
            // the following is copied from TTree::SetBranchAddress<T>
            b.ROOTClass = TClass::GetClass(typeid(T));
            b.ROOTType = kOther_t;
            if (b.ROOTClass==0) b.ROOTType = TDataType::GetType(typeid(T));
            return b;
        }
        std::string Name;
        TClass* ROOTClass;
        EDataType ROOTType;
        void** ValuePtr;
    };

    std::vector<ROOT_branch_t> branches;
};

}

// macro to define branches consistently
#define ADD_BRANCH_T(type, name) Branch_t<type> name{this, #name};
