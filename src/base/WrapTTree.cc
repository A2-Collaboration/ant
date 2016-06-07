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
    auto tree_trick = (TTree_trick*)Tree;
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
    vector<string> branch_names;
    for(int i=0;i<tree->GetNbranches();i++) {
        TBranch* b = dynamic_cast<TBranch*>(tree->GetListOfBranches()->At(i));
        branch_names.emplace_back(b->GetName());
    }

    // ensure exact match of branches if required
    if(exact && branch_names.size() != branches.size()) {
        LOG_IF(!nowarn, WARNING) << "TTree has wrong number of branches";
        return false;
    }

    // ensure that we find all branches
    for(const ROOT_branch_t& b : branches) {
        /// \todo one could check types here as well, not only names
        if(!std_ext::contains(branch_names, b.Name)) {
            LOG(WARNING) << "Branch " << b.Name << " not found";
            return false;
        }
    }

    // all branch names found, so we match
    return true;
}
