#include "TreeManager.h"
#include "TTree.h"

using namespace ant::analysis::input;

TreeManager::TreeManager()
{
}

TreeManager::~TreeManager()
{
}

void TreeManager::AddTree(TTree *tree)
{
    if(!trees.empty()) {
        if(tree->GetEntries() != GetEntries())
            throw false;
    }

    trees.insert(tree);
}

Long64_t TreeManager::GetEntries() const
{
    const auto& it = trees.begin();
    return it != trees.end() ? (*it)->GetEntries() : 0 ;

}

void TreeManager::GetEntry(Long64_t entry)
{
    for(auto& tree : trees) {
        tree->GetEntry(entry);
    }
}
