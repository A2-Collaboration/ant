#include "WrapTTree.h"

#include "base/std_ext/vector.h"
#include "base/std_ext/string.h"
#include "base/std_ext/memory.h"
#include "base/Logger.h"

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

struct WrapTTree::ROOTArrayNotifier_t : TObject {
    std::list<std::function<void()>> LoadNotifiers;
    virtual bool Notify() override {
        // notify all subscribers
        for(auto& n : LoadNotifiers)
            n();
        // TTree/TChain don't check this return value...
        // but indicate success wrap anyway...
        return true;
    }
};

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
        // the default, some very simple type
        return t.SetBranchAddress(b.Name.c_str(),*b.ValuePtr,0,b.ROOTClass,b.ROOTType,false);
    };

    for(const auto& b : branches) {
        if(b.IsROOTArray) {
            HandleROOTArray(b);
            continue;
        }
        const auto res = set_branch_address(*Tree, b);
        if(res < TTree::kMatch)
            throw Exception(std_ext::formatter() << "Cannot set branch " << b.Name << " in tree " << Tree->GetName());
    }

    Tree->SetNotify(ROOTArrayNotifier.get());
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
        EDataType rootType;
        TClass* rootClass;
        if (b->GetExpectedType(rootClass,rootType) != 0) {
            LOG(ERROR) << "Given branch did not tell expected class/type";
            return false;
        }
        treebranches.emplace_back(b->GetName(), rootClass, rootType);
    }

    // ensure exact match of branches if required
    if(exact && treebranches.size() != branches.size()) {
        LOG_IF(!nowarn, WARNING) << "TTree has wrong number of branches";
        return false;
    }

    // ensure that we find all branches
    for(const ROOT_branch_t& b : branches) {
        auto it_treebranch = std::find(treebranches.begin(), treebranches.end(), b.Name);
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
        auto it_b = std::find(branches.begin(), branches.end(), src_b.Name);
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

WrapTTree::WrapTTree() :
    ROOTArrayNotifier(std_ext::make_unique<ROOTArrayNotifier_t>())
{}

WrapTTree::~WrapTTree()
{
    // make we don't linger around in TTree fNotify after destruction
    // NOTE: If you encounter segfaults here, the TTree was destroyed
    // (for example because owning TFile closed) before this WrapTTree
    // was destroyed.
    // Fix: Change order of initialization/destruction of WrapTFileInput and WrapTTree
    // (see UnpackerA2Geant as example, commit d7efcda8f)
    if(Tree && Tree->GetNotify() == ROOTArrayNotifier.get())
        Tree->SetNotify(0);
}

void WrapTTree::HandleROOTArray(const WrapTTree::ROOT_branch_t& b)
{
    // handle this quite specially
    auto ROOT_branch = Tree->GetBranch(b.Name.c_str());

    if(!ROOT_branch)
        throw Exception("WrapTTree::Array: Cannot find  branch "+b.Name+" in tree");
    // presumably splitted branches have more than one leaf?
    if(ROOT_branch->GetNleaves() != 1)
        throw Exception("WrapTTree::Array: Branch "+b.Name+" does not have exactly one leaf");

    struct TLeaf_wrapper : TLeaf {
        TLeaf_wrapper(TLeaf* leaf, ROOTArray_traits& array) :
            TLeaf(*leaf), // call the copy ctor to setup this wrapper properly
            Leaf(leaf),
            Array(array)
        {
            // "fixed" size branches are easier to handle
            if(!Leaf->GetLeafCount()) {
                if(fLen<=0)
                    throw Exception("Found ROOTArray with no LeafCount branch and length<=0");
                Array.ROOTArray_setSize(fLen);
                Leaf->SetAddress(Array.ROOTArray_getPtr());
            }
            else {
                if(fLen != 1)
                    throw Exception("Found multi-dim ROOTArray, not supported yet");
            }
        }

        TLeaf* const Leaf;
        ROOTArray_traits& Array;

        virtual void ReadBasket(TBuffer& b) override {
            // the code is inspired TLeafD::ReadBasket

            if(fLeafCount) {
                // make sure the leafCount branch looks at the same entry as this branch
                // we do this to know the size beforehand
                const Long64_t entry = fBranch->GetReadEntry();
                if (fLeafCount->GetBranch()->GetReadEntry() != entry) {
                    fLeafCount->GetBranch()->GetEntry(entry);
                }
                // oh well, what a nice cast here...
                Array.ROOTArray_setSize(Int_t(fLeafCount->GetValue()));

                // (re)set the pointer to our Array, then Leaf's ReadBasket
                // handles the rest (such as byte ordering, arg...)
                // this reset must happen before every read, as ROOTArray_setSize
                // might have cause re-allocation of Array
                Leaf->SetAddress(Array.ROOTArray_getPtr());
            }

            Leaf->ReadBasket(b);
        }

        virtual ~TLeaf_wrapper() {
            // properly clean-up our adopted child
            delete Leaf;
        }
    };

    // as there's only one leaf, it's easy to replace by wrapper
    // however, we need to do this via a "load notifier", as we might deal with TChain
    // here which constantly creates/destroys TTrees (and associated branches and leafs)
    // Notify is called after the TTree (of a TChain, maybe) was loaded

    struct LoadNotifier : TObject {
        TTree& Tree;
        ROOTArray_traits& Array;
        const string BranchName;

        LoadNotifier(TTree& tree, ROOTArray_traits& array, const string& branchName)
            : Tree(tree), Array(array), BranchName(branchName)
        {
            // execute Notify at least once, since if we're not a TChain,
            // TTree might already be loaded when calling LinkBranches()
            // (calling GetEntries() or so is enough to trigger TTree::Load)
            operator()();
        }

        // make it callable
        void operator()() const {
            auto ROOT_branch = Tree.GetBranch(BranchName.c_str());
            auto ROOT_leaf = ROOT_branch->GetListOfLeaves()->GetObjectRef();
            // check if this leaf was already wrapped
            // (appears to be never run by test, but better be safe than sorry)
            if(dynamic_cast<TLeaf_wrapper*>(*ROOT_leaf) != nullptr) {
                return;
            }
            *ROOT_leaf= new TLeaf_wrapper(
                            dynamic_cast<TLeaf*>(*ROOT_leaf),
                            Array
                            );
        }
    };

    ROOTArrayNotifier->LoadNotifiers.emplace_back(LoadNotifier(
                *Tree,
                // ugly cast to get our Array traits back
                *static_cast<ROOTArray_traits*>(*b.ValuePtr),
                b.Name
                ));
}
