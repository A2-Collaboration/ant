#ifndef TREEMANAGER_H
#define TREEMANAGER_H

#include <list>
#include "Rtypes.h"

class TTree;

namespace ant {
namespace input {

class TreeManager {
protected:
    using tree_list_t = std::list<TTree*>;
    tree_list_t trees;

public:
    TreeManager();
    virtual ~TreeManager();

    virtual void AddTree(TTree* tree);

    virtual Long64_t GetEntries() const;

    virtual void GetEntry(Long64_t entry);

};
}
}

#endif
