#include "tree/TEvent.h"
#include "TFile.h"
#include "TTree.h"
#include "analysis/input/ant/detail/Convert.h"
#include "analysis/data/Event.h"
#include <iostream>

using namespace std;
using namespace ant;


int main() {

    TFile f2("TEvent_test.root","READ");

    if(!f2.IsOpen()) {
        cerr << "can't ope file" << endl;
        return 1;
    }

    TTree* tree(nullptr);

    f2.GetObject("teventtest", tree);

    if(tree==nullptr) {
        cerr << "Tree not found" << endl;
        return 2;
    }

    TEvent* readback = nullptr;

    tree->SetBranchAddress("event",&readback);

    if(tree->GetEntries()!=1)
        return 2;

    tree->GetEntry(0);

    cout << *readback << endl;

    auto antevent = ant::input::Convert(*readback);

    cout << *antevent << endl;

    return 0;
}


