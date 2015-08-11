/**
  * @file test_display.cc
  * @breif Demo for Setmarker() of the cbtaps display. can be deleted soon
  * @note crashes on exit. no idea.
  */

#include "TCanvas.h"
#include "TRint.h"
#include "base/cbtaps_display/TH2CB.h"
#include "base/cbtaps_display/TH2TAPS.h"
#include "TMarker.h"

using namespace ant;

int main(int argc, char** argv) {
    TRint app("",&argc,argv);

    TH2CB* h = new TH2CB();
    h->FillElementNumbers();
    h->Draw("col");

    auto m = h->SetMarker(0);
    m->SetMarkerColor(kWhite);
    m->Draw();

    auto m2 = h->SetMarker(1);
    m2->SetMarkerColor(kPink);
    m2->Draw();

    auto m3 = h->SetMarker(719);
    m3->SetMarkerColor(kBlue);
    m3->Draw();

    app.Run();
}
