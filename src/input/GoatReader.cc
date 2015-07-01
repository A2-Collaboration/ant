#include "GoatReader.h"

#include <string>
#include <iostream>

#include "TTree.h"

using namespace ant;
using namespace ant::input;
using namespace std;

void GoatReader::AddInputFile(const std::string &filename)
{
    files.OpenFile(filename);
}

class MyTreeRequestMgr: public TreeRequestManager {
protected:
    FileManager& m;
    TreeManager& t;

public:
    MyTreeRequestMgr(FileManager& _m, TreeManager& _t):
        m(_m), t(_t) {}

    TTree *GetTree(const std::string &name) {
        TTree* tree = nullptr;
        if( m.GetObject(name, tree) ) {
            t.AddTree(tree);
            return tree;
        } else
            return nullptr;
    }

};

//TODO: find a smart way to manage trees and modules:
//   if module does not init or gets removed-> remove also the tree from the list
//   two modules use same tree?
//   reset branch addresses ?

void GoatReader::Initialize()
{
    for(auto module = active_modules.begin(); module != active_modules.end(); ) {

        if( (*module)->SetupBranches( MyTreeRequestMgr(files, trees))) {
            module++;
        } else {
            module = active_modules.erase(module);
        }
    }

}

Long64_t GoatReader::GetNEvents() const
{
    return trees.GetEntries();
}

bool GoatReader::hasData() const {
    return current_entry < GetNEvents();
}

std::shared_ptr<Event> GoatReader::ReadNextEvent()
{
    ++current_entry;
    trees.GetEntry(current_entry);

    active_modules.GetEntry();

    cout << "Entry " << current_entry << endl;
    cout << "Tagger hits: " << tagger.GetNTagged() << endl;

    std::shared_ptr<Event> event = make_shared<Event>();

    return event;
}
