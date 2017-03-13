#include "WrapTTree.h"

#include "base/std_ext/vector.h"
#include "base/Logger.h"
#include "base/std_ext/string.h"

#include "TBufferFile.h"
#include "TLeaf.h"

using namespace std;
using namespace ant;

void WrapTTree::CreateBranches(TTree* tree) {
    // some checks first
    if(tree==nullptr)
        throw Exception("Provided null tree to CreateBranches");

    for(const auto& b : branches) {
        if(b.IsROOTArray)
            throw Exception("Cannot use ROOTArray branch "+b.Name+" for output");
    }

    Tree = tree;
    // little trick to access the protected method
    struct TTree_trick : TTree {
        using TTree::BranchImpRef;
    };
    auto tree_trick = reinterpret_cast<TTree_trick*>(Tree);
    for(const auto& b : branches) {
        // dereference ValuePtr here to pointer to value
        tree_trick->BranchImpRef(b.Name.c_str(), b.ROOTClass, b.ROOTType, *b.ValuePtr, 32000, 99);
    }
}

int WrapTTree::ROOT_branch_t::HandleROOTArray(TTree& t) const
{
    // handle this quite specially...
    auto ROOT_branch = t.GetBranch(Name.c_str());

    if(!ROOT_branch)
        throw Exception("WrapTTree::Array: Cannot find  branch "+Name+" in tree");
    if(ROOT_branch->GetNleaves() != 1)
        throw Exception("WrapTTree::Array: Splitted branches are not supported (nLeaves != 1)");

    struct TLeaf_wrapper : TLeaf {
        TLeaf_wrapper(TLeaf* leaf, ROOTArray_traits& array) :
            TLeaf(*leaf), // call the copy ctor to setup this wrapper properly
            Leaf(leaf),
            Array(array)
        {
            // require that this branch is of type "name[size]/D"
            // when it was created with TTree::Branch
            // other cases are not handled by this
            if(!Leaf->GetLeafCount()) {
                throw Exception("The leaf "+string(GetName())+" does not have an associated LeafCount");
            }
        }

        TLeaf* const Leaf;
        ROOTArray_traits& Array;

        virtual void ReadBasket(TBuffer& b) override {
            // assumes fLeafCount is not null, see ctor
            // the code is based on TLeafD::ReadBasket

            // make sure the leafCount branch looks at the same entry as this branch
            // we do this to know the size beforehand
            Long64_t entry = fBranch->GetReadEntry();
            if (fLeafCount->GetBranch()->GetReadEntry() != entry) {
               fLeafCount->GetBranch()->GetEntry(entry);
            }
            // oh well, what a nice cast here...
            Int_t len = Int_t(fLeafCount->GetValue());

            Array.ROOTArray_setSize(len*fLen);

            // set the pointer to our Array, then Leaf's ReadBasket
            // handles the rest (such as byte ordering, arg...)
            Leaf->SetAddress(Array.ROOTArray_getPtr());

            Leaf->ReadBasket(b);
        }

        virtual ~TLeaf_wrapper() {
            // properly clean-up our adopted child
            delete Leaf;
        }
    };

    // as there's only one leaf, it's easy to replace by wrapper
    auto ROOT_leaf = ROOT_branch->GetListOfLeaves()->GetObjectRef();
    *ROOT_leaf= new TLeaf_wrapper(
                    // it's not ROOT if you don't use ugly casts...sigh
                    dynamic_cast<TLeaf*>(*ROOT_leaf),
                    *static_cast<ROOTArray_traits*>(*ValuePtr)
                    );

    // successful match
    return static_cast<int>(TTree::kMatch);
}

void WrapTTree::LinkBranches(TTree* tree) {
    if(tree != nullptr)
        Tree = tree;
    if(!Tree)
        throw Exception("Set the Tree pointer (or provide as argument) before calling LinkBranches");

    auto set_branch_address = [] (TTree& t, const ROOT_branch_t& b) {
        // logic copied from TTree::SetBranchAddress<T>
        if(b.ROOTClass) {
            return t.SetBranchAddress(b.Name.c_str(),b.ValuePtr,0,b.ROOTClass,b.ROOTType,true);
        }
        if(b.IsROOTArray) {
            return b.HandleROOTArray(t);
        }
        // the default, some very simple type
        return t.SetBranchAddress(b.Name.c_str(),*b.ValuePtr,0,b.ROOTClass,b.ROOTType,false);
    };

    for(const auto& b : branches) {
        const auto res = set_branch_address(*Tree, b);
        if(res < TTree::kMatch)
            throw Exception(std_ext::formatter() << "Cannot set branch " << b.Name << " in tree " << Tree->GetName());
    }
}

bool WrapTTree::Matches(TTree* tree, bool exact, bool nowarn) const {
    if(tree == nullptr)
        tree = Tree;
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

bool WrapTTree::CopyFrom(const WrapTTree& src) {
    // src has more branches than us, this can't work
    if(src.branches.size()>branches.size())
        return false;

    // copy branches by name
    for(const ROOT_branch_t& src_b : src.branches) {
        auto it_b = std::find(branches.begin(), branches.end(), src_b);
        // src branch not found in our list of branches
        if(it_b == branches.end())
            return false;
        if(it_b->ROOTType != src_b.ROOTType)
            return false;
        // for complex types, we rely on the ROOT machinery to copy them
        if(src_b.ROOTType == kOther_t) {
            if(it_b->ROOTClass != src_b.ROOTClass)
                return false;
            TBufferFile R__b(TBufferFile::EMode::kWrite);
            R__b.WriteClassBuffer(src_b.ROOTClass, *src_b.ValuePtr);
            R__b.SetBufferOffset(); // rewind buffer to read from it
            R__b.ReadClassBuffer(it_b->ROOTClass,*(it_b->ValuePtr),0);
        }
        else {
            auto datatype = TDataType::GetDataType(src_b.ROOTType);
            std::memcpy(*(it_b->ValuePtr), *src_b.ValuePtr, datatype->Size());
        }
    }

    return true;
}
