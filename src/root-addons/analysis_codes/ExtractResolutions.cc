#include "ExtractResolutions.h"

#include "TreeTools.h"
#include "base/std_ext/string.h"
#include "base/interval.h"

#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "base/std_ext/math.h"
#include "TF1.h"
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::std_ext;

double LogNormalPDF(const double* X, const double* p) {

    const auto& x  = X[0];

    const auto& A  = p[0];
    const auto& mu = p[1];
    const auto& sigma = p[2];

    return A / ( x * sigma * sqrt(2*M_PI) ) * exp( - sqr(log(x) - mu) / (2 * sqr(sigma)));
}

void ExtractResolutions::AnalyseThetaCB(TTree* tree)
{
    const interval<double> range = {-.2,.2};
    const unsigned nElements = 720;

    new TCanvas();
    auto h = Draw(tree, "Theta:Element",TCut(""), int(nElements), 0, nElements, 80, range.Start(), range.Stop());
    h->SetXTitle("Element");
    h->SetYTitle("#Delta #theta [rad]");
    h->SetZTitle("# Particles");

    TH2CB* cb = new TH2CB("theta_sigmas", "theta sigmas");
    cb->SetZTitle("#sigma_{#theta} [rad]");

    TF1* f = new TF1("f", "gaus", range.Start(), range.Stop());

    new TCanvas();
    for(unsigned i=0; i<nElements; ++i) {

        f->SetParameter(0, 1);
        f->SetParameter(1, 0);
        f->SetParameter(2,.1);
        TH1* pr;
        pr = h->ProjectionY(Form("Theta_Element_%d",i), i, i+1);
        pr->SetTitle(Form("Element %d", i));
        pr->SetXTitle("#Delta #theta [rad]");
        pr->Draw();
        pr->Fit(f, "RQ");
        cb->SetElement(i, f->GetParameter(2));

//        gPad->SaveAs(Form("Theta_cb_%d.png",i));
//        gPad->SaveAs(Form("Theta_cb_%d.pdf",i));
    }

    new TCanvas();
    cb->Draw("colz");
}

void ExtractResolutions::AnalysePhiCB(TTree* tree)
{
    const interval<double> range = {-.2,.2};
    const unsigned nElements = 720;

    new TCanvas();
    auto h = Draw(tree, "Phi:Element",TCut(""), int(nElements), 0, nElements, 80, range.Start(), range.Stop());
    h->SetXTitle("Element");
    h->SetYTitle("#Delta #phi [rad]");
    h->SetZTitle("# Particles");

    TH2CB* cb = new TH2CB("phi_sigmas", "phi sigmas");
    cb->SetZTitle("#sigma_{#phi} [rad]");

    TF1* f = new TF1("f", "gaus", range.Start(), range.Stop());

    new TCanvas();
    for(unsigned i=0; i<nElements; ++i) {

        f->SetParameter(0, 1);
        f->SetParameter(1, 0);
        f->SetParameter(2,.1);
        TH1* pr;
        pr = h->ProjectionY(Form("Phi_Element_%d",i), i, i+1);
        pr->SetTitle(Form("Element %d", i));
        pr->SetXTitle("#Delta #phi [rad]");
        pr->Draw();
        pr->Fit(f, "RQ");
        cb->SetElement(i, f->GetParameter(2));

//        gPad->SaveAs(Form("Phi_cb_%d.png",i));
//        gPad->SaveAs(Form("Phi_cb_%d.pdf",i));
    }

    new TCanvas();
    cb->Draw("colz");
}

void ExtractResolutions::AnalyseECB(TTree* tree)
{
    const interval<double> range = {-120,0};
    const unsigned nElements = 720;

    new TCanvas();
    auto h = Draw(tree, "E:Element",TCut(""), int(nElements), 0, nElements, 120, range.Start(), range.Stop());
    h->SetXTitle("Element");
    h->SetYTitle("#Delta E [MeV]");
    h->SetZTitle("# Particles");

    TH2CB* cb = new TH2CB("E_sigmas", "E sigmas");
    cb->SetZTitle("#sigma_{E} [MeV]");

    TF1* f = new TF1("f", "gaus", range.Start(), range.Stop());

    new TCanvas();
    for(unsigned i=0; i<nElements; ++i) {

        f->SetParameter(0, 1);
        f->SetParameter(1, 0);
        f->SetParameter(2,.1);
        TH1* pr;
        pr = h->ProjectionY(Form("E_Element_%d",i), i, i+1);
        pr->SetTitle(Form("Element %d", i));
        pr->SetXTitle("#Delta E [MeV]");
        pr->Draw();
        pr->Fit(f, "RQ");
        cb->SetElement(i, f->GetParameter(2));

//        gPad->SaveAs(Form("Phi_cb_%d.png",i));
//        gPad->SaveAs(Form("Phi_cb_%d.pdf",i));
    }

    new TCanvas();
    cb->Draw("colz");
}

TF1*ExtractResolutions::LogNormal()
{
    return new TF1("lognormal", LogNormalPDF, 0, 10, 3);
}
