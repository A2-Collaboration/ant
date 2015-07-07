#include "input/GoatReader.h"
#include "data/Event.h"

#include <iostream>

using namespace std;
using namespace ant;

int main(int argc, char** argv) {

    ant::input::GoatReader g;

    for(int i = 1; i < argc; ++i) {
        g.AddInputFile(argv[i]);
    }

    g.Initialize();

    unsigned int n = 0;
    while(g.hasData() && (n++ < 10)) {
        auto event = g.ReadNextEvent();
        cout << *event << endl;
    }

    return 0;
}
