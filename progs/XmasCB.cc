#include "base/Logger.h"
#include "TRint.h"
#include "base/CmdLine.h"
#include "base/cbtaps_display/TH2CB.h"
#include "TCanvas.h"

using namespace std;
using namespace ant;

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("XmasCB", ' ', "0.1");
    cmd.parse(argc, argv);

    int fake_argc=0;
    char** fake_argv=nullptr;
    TRint app("Ant",&fake_argc,fake_argv,nullptr,0,true);

    const int w_px = 1024;

    TH2CB* cb = new TH2CB("cb","");
    cb->SetTitle("");

    TCanvas* c = new TCanvas("ccb","CB");
    c->SetMargin(0,0,0,0);

    const double x1 = cb->GetXaxis()->GetXmin();
    const double x2 = cb->GetXaxis()->GetXmax();
    const double y1 = cb->GetYaxis()->GetXmin();
    const double y2 = cb->GetYaxis()->GetXmax();

    const double ratio = (x2-x1)/(y2-y1);

    c->SetCanvasSize(w_px, unsigned(w_px/ratio));
    c->cd();
    cb->Draw("col");

    TH2CB* grid = new TH2CB("grid");
    grid->SetTitle("");
    grid->Draw("same");


    app.Run();


    return 0;
}
