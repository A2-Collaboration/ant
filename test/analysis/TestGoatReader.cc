#include "catch.hpp"

#include "analysis/input/goat/GoatReader.h"
#include "analysis/data/Event.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();


TEST_CASE("GoatReader", "[analysis]") {
    dotest();
}

void dotest() {

    ant::input::GoatReader g;

    /// \todo Generate or read some Goat file?

    g.AddInputFile("NOTTHEREYET");

    g.Initialize();

    unsigned int n = 0;
    while(g.hasData() && (n++ < 10)) {
        auto event = g.ReadNextEvent();
        cout << *event << endl;
    }
}
