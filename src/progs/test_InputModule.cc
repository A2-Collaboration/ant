#include "input/FileManager.h"
#include "input/TreeManager.h"
#include "input/PlutoInput.h"
#include "TTree.h"

#include <iostream>

using namespace ant::input;
using namespace std;

class MyTreeRequestMgr: public TreeRequestManager {
protected:
    FileManager& m;
    TreeManager& t;

public:
    MyTreeRequestMgr(FileManager& _m, TreeManager& _t):
        m(_m), t(_t) {}

    TTree *GetTree(const string &name) {
        TTree* tree = nullptr;
        if( m.GetObject(name, tree) ) {
            t.AddTree(tree);
            return tree;
        } else
            return nullptr;
    }
};

int main(int argc, char** argv) {

    FileManager m;

    for(int i = 1; i < argc; ++i) {
        m.OpenFile(argv[i]);
    }

    TreeManager treeManager;

    PlutoInput pluto;

    bool init = pluto.SetupBranches( MyTreeRequestMgr(m,treeManager) );

    if( init == true) {
        cout << "Init OK" << endl;

        cout << treeManager.GetEntries() << " entries " << endl;

        treeManager.GetEntry(1);
        pluto.GetEntry();

        cout << "Pluto particles: " << pluto.Particles().size() << endl;
    }

    m.CloseAll();

}
