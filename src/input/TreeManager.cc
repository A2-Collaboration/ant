#include "TreeManager.h"
#include "TTree.h"

using namespace ant;
using namespace input;

TreeManager::TreeManager()
{

}

ant::input::TreeManager::~TreeManager()
{

}

void ant::input::TreeManager::AddTree(TTree *tree)
{
    if(!trees.empty()) {
        if(tree->GetEntries() != GetEntries())
            throw false;
    }

    trees.push_back(tree);
}

Long64_t ant::input::TreeManager::GetEntries() const
{
    return trees.front()->GetEntries();
}

void ant::input::TreeManager::GetEntry(Long64_t entry)
{
    for(auto& tree : trees) {
        tree->GetEntry(entry);
    }
}
