#include "calibration/gui/ACECanvas.h"
#include "calibration/gui/Dialogs.h"
#include <iostream>
#include "TRint.h"

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;


int main(int argc, char** argv) {
//    SetupLogger();

    int i = 0;
    TRint app("ACE",&i,nullptr);

    if (argc == 2)
        new ACECanvas(argv[1]);
    else
    {
        cerr << "Usage: " << argv[0] << " <calibration-file>" << endl;
        return 1;
    }

    app.Run(kFALSE);

}
