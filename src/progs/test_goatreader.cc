#include "input/GoatReader.h"

#include <iostream>

int main(int argc, char** argv) {

    ant::input::GoatReader g;

    for(int i = 1; i < argc; ++i) {
        g.AddInputFile(argv[i]);
    }

    g.Initialize();

    while(g.hasData()) {
        g.ReadNextEvent();
    }

    return 0;
}
