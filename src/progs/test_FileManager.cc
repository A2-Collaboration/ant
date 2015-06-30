#include "input/FileManager.h"
#include "TTree.h"

#include <iostream>

using namespace ant::input;
using namespace std;

int main(int argc, char** argv) {

    FileManager m;

    for(int i = 1; i < argc; ++i) {
        m.OpenFile(argv[i]);
    }

    TTree* tree = nullptr;

    if( m.GetObject("data", tree) ) {
        cout << "Found tree " << tree->GetName() << endl;
        cout << tree->GetEntries() << endl;
    }

    m.CloseAll();

}
