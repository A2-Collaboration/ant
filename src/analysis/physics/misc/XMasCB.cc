#include "XMasCB.h"
#include "base/cbtaps_display/TH2CB.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace ant;
using namespace ant::analysis::physics;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace std;

XMasCB::XMasCB(const std::string& name, PhysOptPtr opts):
    Physics(name, opts),
    ext(opts->Get<string>("FileType", "pdf")),
    w_px(opts->Get<int>("x_px", 1024))
{
    hist = std_ext::make_unique<TH2CB>("", "", opts->Get<bool>("GluePads", true));
    hist->SetTitle("");
    c    = new TCanvas("xmascb", "XMasCB");
}

XMasCB::~XMasCB()
{}

void XMasCB::ProcessEvent(const Event& event)
{
    hist->ResetElements(0.0);

    for(const auto& c : event.Reconstructed().AllClusters()) {
        if(c.Detector & Detector_t::Type_t::CB) {
            for(const auto& hit : c.Hits) {
                for(const auto& datum : hit.Data) {
                    if(datum.Type == Channel_t::Type_t::Integral) {
                        hist->SetElement(hit.Channel, hist->GetElement(hit.Channel)+datum.Value);
                    }
                }
            }
        }
    }

    c->cd();
    c->SetLogz();
    c->SetMargin(0,0,0,0);
    //gStyle->SetPalette(1);
    //hist->SetContour(n);

    // comment: this color palette has three colors, therefore the arrays for red,
    // green and blue have three entries. The first entry in *every* array combined
    // describes the first color of the palette. So if you want red as first color,
    // the first entry of red has to be 1, of green and blue 0. Rather unintuitive,
    // but at least it works.
    const Int_t Number = 3;

    //gelb-grün-blau
//    Double_t Red[Number]    = { 1, 0.27, 0.0};
//    Double_t Green[Number]  = { 0.82, 0.86, 0.44};
//    Double_t Blue[Number]   = { 0.34, 0.26, 0.87};
    //gelb-rot-violett
//    Double_t Red[Number]    = { 1, 0.77, 0.53};
//    Double_t Green[Number]  = { 0.8, 0.16, 0.03};
//    Double_t Blue[Number]   = { 0.23, 0.06, 0.91};
    //hellblau-magenta-orange
    Double_t Red[Number]    = { 0.56, 0.9, 1};
    Double_t Green[Number]  = { 0.66, 0.02, 0.65};
    Double_t Blue[Number]   = { 1, 0.44, 0};

    Double_t Length[Number] = { 0.00, 0.10, 1.00 };
    Int_t nb=50;
    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
    hist->SetContour(nb);

    const double x1 = hist->GetXaxis()->GetXmin();
    const double x2 = hist->GetXaxis()->GetXmax();
    const double y1 = hist->GetYaxis()->GetXmin();
    const double y2 = hist->GetYaxis()->GetXmax();

    const double ratio = (x2-x1)/(y2-y1);

    c->SetCanvasSize(w_px, unsigned(w_px/ratio));

    hist->Draw("col");

    std_ext::formatter f;
    f << "xmas_cb_" << n++ << "." << ext;

    c->SaveAs(f.str().c_str());

}

AUTO_REGISTER_PHYSICS(XMasCB)
