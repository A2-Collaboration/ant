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

    trees.insert(tree);
}

Long64_t ant::input::TreeManager::GetEntries() const
{
    const auto& it = trees.begin();
    return it != trees.end() ? (*it)->GetEntries() : 0 ;

}

void ant::input::TreeManager::GetEntry(Long64_t entry)
{
    for(auto& tree : trees) {
        tree->GetEntry(entry);
    }
}
