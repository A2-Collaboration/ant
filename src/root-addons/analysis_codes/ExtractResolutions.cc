#include "ExtractResolutions.h"

#include "TreeTools.h"
#include "analysis/plot/Histogram.h"
#include "base/std_ext/string.h"
#include "base/interval.h"

#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TDirectory.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "root-addons/cbtaps_display/TH2TAPS.h"
#include "base/std_ext/math.h"
#include "TF1.h"
#include <cmath>
#include "analysis/plot/root_draw.h"

#include "analysis/utils/KinFitter.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;

void Analyse(TTree* tree, const unsigned nElements, const string& branch, const string& detname, const BinSettings energy_bins, const BinSettings vbins);

double LogNormalPDF(const double* X, const double* p) {

    const auto& x  = X[0];

    const auto& A  = p[0];
    const auto& mu = p[1];
    const auto& sigma = p[2];

    return A / ( x * sigma * sqrt(2*M_PI) ) * exp( - sqr(log(x) - mu) / (2 * sqr(sigma)));
}

void ExtractResolutions::AnalyseThetaCB(TTree* tree)
{
    Analyse(tree, 720, "Phi", "CB", BinSettings(16,0,1600), BinSettings(80,-.2,.2));
}

void ExtractResolutions::AnalysePhiCB(TTree* tree)
{
    Analyse(tree, 720, "Phi", "CB", BinSettings(16,0,1600), BinSettings(80,-.2,.2));
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
        pr = h->ProjectionY(Form("E_Element_%d",i), int(i), int(i+1));
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

TF1*ExtractResolutions::SigmaFit()
{
    return utils::KinFitter::angular_sigma::GetTF1();
}

void ExtractResolutions::AnalyseThetaTAPS(TTree* tree) {
    Analyse(tree, 438, "Theta", "TAPS", BinSettings(16,0,1600), BinSettings(80,-.2,.2));
}

void ExtractResolutions::AnalysePhiTAPS(TTree* tree) {
    Analyse(tree, 438, "Phi", "TAPS", BinSettings(16,0,1600), BinSettings(80,-.5,.5));
}

void Analyse(TTree* tree, const unsigned nElements, const string& branch, const string& detname, const BinSettings energy_bins, const BinSettings vbins)
{
    const BinSettings xbins(nElements, 0, nElements);

    TCanvas* c = new TCanvas();
    auto h = Draw(tree, Form("%s:tE:Element", branch.c_str()), TCut(""), xbins, energy_bins, vbins);
    h->SetXTitle("Element");
    h->SetZTitle(Form("#Delta %s [rad]",branch.c_str()));
    h->SetYTitle("E [MeV]");
    h->SetName(Form("%s_%s", detname.c_str(), branch.c_str()));

    TF1* f = new TF1("f", "gaus", vbins.Start(), vbins.Stop());
    f->SetParameter(0, 1);
    f->SetParameter(1, 0);
    f->SetParameter(2,.1);



    TCanvas* fits = new TCanvas();
    bool first = true;
    c->cd();


    TH1D* h_global_sigma = new TH1D(Form("%s_%s_global_sigma",    detname.c_str(), branch.c_str()), "Global #Delta #theta, taps",  int(nElements), 0, nElements);
    TH1D* h_fit_p0       = new TH1D(Form("%s_%s_global_sigma_p0", detname.c_str(), branch.c_str()), "Fit #sigma #theta, p0, taps", int(nElements), 0, nElements);
    TH1D* h_fit_p1       = new TH1D(Form("%s_%s_global_sigma_p1", detname.c_str(), branch.c_str()), "Fit #sigma #theta, p1, taps", int(nElements), 0, nElements);
    TH1D* h_fit_p2       = new TH1D(Form("%s_%s_global_sigma_p2", detname.c_str(), branch.c_str()), "Fit #sigma #theta, p2, taps", int(nElements), 0, nElements);

    for(int e=0; e<int(nElements); ++e) {
        h->GetXaxis()->SetRange(e,e+1);

        const string projcmd = to_string(e)+"_zy";

        TH2* proj = dynamic_cast<TH2*>(h->Project3D(projcmd.c_str()));
        proj->SetTitle(Form("#Delta %s, %s, Element %d", branch.c_str(), detname.c_str(),e));
        proj->Draw("colz");

        proj->FitSlicesY(f,0,-1,0,"QNR");

        TH1D* h_sigma = nullptr;

        gDirectory->GetObject(Form("%s_%s_%d_zy_2", detname.c_str(),branch.c_str(),e), h_sigma);
        h_sigma->SetTitle(Form("#Delta %s, fitted, %s, Element %d",branch.c_str(), detname.c_str(),e));

        TF1* fe = utils::KinFitter::angular_sigma::GetTF1(Form("f_%s_%s_%d",detname.c_str(),branch.c_str(), e));
        h_sigma->Fit(fe,"MRQ");

        TH1* p1 = proj->ProjectionY(Form("%s_%s_%d_proj",detname.c_str(),branch.c_str(),e),0,-1);
        TF1* p1f = new TF1("", "gaus", vbins.Start(), vbins.Stop());
        p1f->SetParameter(0, 1);
        p1f->SetParameter(1, 0);
        p1f->SetParameter(2,.1);
        p1->SetTitle(Form("#Delta %s, global, %s, Element %d",branch.c_str(), detname.c_str(),e));
        p1->Fit(p1f,"QNR");

        cout << "Element " << e << ":  Global sigma="<< p1f->GetParameter(2) <<", E-Fitted: "
             << fe->GetParameter(0) << " "
             << fe->GetParameter(1) << " "
             << fe->GetParameter(2)
             << " Chi2/dof=" << fe->GetChisquare()/fe->GetNDF() << endl;

        h_global_sigma->SetBinContent(e+1,p1f->GetParameter(2));
        h_fit_p0->SetBinContent(e+1, fe->GetParameter(0));
        h_fit_p1->SetBinContent(e+1, fe->GetParameter(1));
        h_fit_p2->SetBinContent(e+1, fe->GetParameter(2));

        fits->cd();
        if(first) {
            fe->Draw();
            first=false;
        } else {
            fe->Draw("same");
        }
        c->cd();


    }


    TH2Crystals* h_global = nullptr;
    TH2Crystals* h_p0 = nullptr;
    TH2Crystals* h_p1 = nullptr;
    TH2Crystals* h_p2 = nullptr;

    if(detname == "CB") {
        h_global = new TH2CB(Form("%s_sigmas",    branch.c_str()), Form("#sigma %s, global, %s", branch.c_str(), detname.c_str()));
        h_p0     = new TH2CB(Form("%s_sigmas_p0", branch.c_str()), Form("%s sigmas p0", branch.c_str()));
        h_p1     = new TH2CB(Form("%s_sigmas_p1", branch.c_str()), Form("%s sigmas p1", branch.c_str()));
        h_p2     = new TH2CB(Form("%s_sigmas_p2", branch.c_str()), Form("%s sigmas p2", branch.c_str()));
    } else if (detname == "TAPS") {
        h_global = new TH2TAPS(Form("%s_sigmas",    branch.c_str()), Form("#sigma %s, global, %s", branch.c_str(), detname.c_str()));
        h_p0     = new TH2TAPS(Form("%s_sigmas_p0", branch.c_str()), Form("%s sigmas p0", branch.c_str()));
        h_p1     = new TH2TAPS(Form("%s_sigmas_p1", branch.c_str()), Form("%s sigmas p1", branch.c_str()));
        h_p2     = new TH2TAPS(Form("%s_sigmas_p2", branch.c_str()), Form("%s sigmas p2", branch.c_str()));
    }

    h_global->SetZTitle(Form("#sigma %s [rad]", branch.c_str()));
    h_global->SetElements(*h_global_sigma);

    h_p0->SetZTitle("p0");
    h_p0->SetElements(*h_fit_p0);

    h_p1->SetZTitle("p1");
    h_p1->SetElements(*h_fit_p1);

    h_p2->SetZTitle("p2");
    h_p2->SetElements(*h_fit_p2);

    canvas("Sigma "+branch+" "+detname)
            << drawoption("colz")
            << h_global
            << h_p0
            << h_p1
            << h_p2
            << endc;

}
