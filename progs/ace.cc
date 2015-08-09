#include "calibration/ACECanvas.h"
#include <iostream>
#include "TRint.h"

using namespace std;
using namespace ant;


//int main(int argc, char** argv) {
int main() {
//    SetupLogger();

    int i = 0;
    TRint app("ACE",&i,nullptr);

    ACECanvas* c = new ACECanvas("test");

    app.Run(kFALSE);

}
