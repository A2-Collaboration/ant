#include "calibration/gui/ManagerWindow.h"

#include "TRint.h"

int main(int argc, char** argv)
{
    TRint app("guitest",&argc,argv);
    // MainFram is destroyed on close!
    new ant::calibration::gui::ManagerWindow(gClient->GetRoot(),400,400);
    app.Run();
}
