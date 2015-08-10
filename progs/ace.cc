#include "calibration/ACECanvas.h"
#include "calibration/ListBox.h"
#include <iostream>
#include "TRint.h"

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;


//int main(int argc, char** argv) {
int main() {
//    SetupLogger();

    int i = 0;
    TRint app("ACE",&i,nullptr);

    new ACECanvas("test.root");
//    new ant::calibration::ListQueryWindow("test","test text",{"a","b","c"});

    app.Run(kFALSE);

}
