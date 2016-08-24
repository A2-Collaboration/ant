#include "WrapTTree.h"

#include "base/std_ext/vector.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;

void WrapTTree::CreateBranches(TTree* tree) {
    Tree = tree;
    // little trick to access the protected method
    struct TTree_trick : TTree {
        using TTree::BranchImpRef;
    };
    auto tree_trick = dynamic_cast<TTree_trick*>(Tree);
    for(const auto& b : branches) {
        // dereference ValuePtr here to pointer to value
        tree_trick->BranchImpRef(b.Name.c_str(), b.ROOTClass, b.ROOTType, *b.ValuePtr, 32000, 99);
    }
}

void WrapTTree::LinkBranches(TTree* tree) {
    Tree = tree;
    for(const auto& b : branches) {
        // copied from TTree::SetBranchAddress<T>
        /// \todo handling of classes should be tested better...
        if(b.ROOTClass)
            Tree->SetBranchAddress(b.Name.c_str(),b.ValuePtr,0,b.ROOTClass,b.ROOTType,true);
        else
            Tree->SetBranchAddress(b.Name.c_str(),*b.ValuePtr,0,b.ROOTClass,b.ROOTType,false);
    }
}

bool WrapTTree::Matches(TTree* tree, bool exact, bool nowarn) const {

    if(!tree)
        return false;
    // get all branches from TTree

    vector<ROOT_branchinfo_t> treebranches;
    for(int i=0;i<tree->GetNbranches();i++) {
        TBranch* b = dynamic_cast<TBranch*>(tree->GetListOfBranches()->At(i));
        treebranches.emplace_back(b->GetName());
        ROOT_branchinfo_t& info = treebranches.back();
        if (b->GetExpectedType(info.ROOTClass,info.ROOTType) != 0) {
            LOG(ERROR) << "Given branch did not tell expected class/type";
            return false;
        }
    }

    // ensure exact match of branches if required
    if(exact && treebranches.size() != branches.size()) {
        LOG_IF(!nowarn, WARNING) << "TTree has wrong number of branches";
        return false;
    }

    // ensure that we find all branches
    for(const ROOT_branch_t& b : branches) {
        auto it_treebranch = std::find(treebranches.begin(), treebranches.end(), b);
        if(it_treebranch == treebranches.end()) {
            LOG_IF(!nowarn, WARNING) << "Did not find branch " << b.Name << " in tree";
            return false;
        }

        if( b.ROOTType != it_treebranch->ROOTType ||
           (b.ROOTType == kOther_t && !b.ROOTClass->InheritsFrom(it_treebranch->ROOTClass) )) {
            LOG_IF(!nowarn, WARNING) << "Branch " << b.Name << " has wrong type";
            return false;
        }
    }

    // all branch names found, so we match
    return true;
}
