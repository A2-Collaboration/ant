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
#include "base/std_ext/math.h"
#include "TF1.h"
#include <cmath>
#include "analysis/plot/root_draw.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;

double LogNormalPDF(const double* X, const double* p) {

    const auto& x  = X[0];

    const auto& A  = p[0];
    const auto& mu = p[1];
    const auto& sigma = p[2];

    return A / ( x * sigma * sqrt(2*M_PI) ) * exp( - sqr(log(x) - mu) / (2 * sqr(sigma)));
}

void ExtractResolutions::AnalyseThetaCB(TTree* tree)
{

    const unsigned nElements = 720;
    const BinSettings xbins(nElements, 0, nElements);
    const BinSettings ybins(16,     0.0, 1600.0);
    const BinSettings zbins(80,    -0.2,    0.2);


    TCanvas* c = new TCanvas();
    auto h = Draw(tree, "Theta:tE:Element", TCut(""), xbins, ybins, zbins);
    h->SetXTitle("Element");
    h->SetZTitle("#Delta #theta [rad]");
    h->SetYTitle("E [MeV]");
    h->SetName("CB_theta");

    TF1* f = new TF1("f", "gaus", zbins.Start(), zbins.Stop());
    f->SetParameter(0, 1);
    f->SetParameter(1, 0);
    f->SetParameter(2,.1);



    TCanvas* fits = new TCanvas();
    bool first = true;
    c->cd();


    TH1D* h_global_sigma = new TH1D("cb_theta_global_sigma", "Global #Delta #theta, cb", nElements, 0, nElements);
    TH1D* h_theta_fit_p0   = new TH1D("cb_theta_fit_p0", "Fit #sigma #theta, p0, cb", nElements, 0, nElements);
    TH1D* h_theta_fit_p1   = new TH1D("cb_theta_fit_p1", "Fit #sigma #theta, p1, cb", nElements, 0, nElements);
    TH1D* h_theta_fit_p2   = new TH1D("cb_theta_fit_p2", "Fit #sigma #theta, p2, cb", nElements, 0, nElements);

    for(int e=0; e<int(nElements); ++e) {
        h->GetXaxis()->SetRange(e,e+1);

        const string projcmd = to_string(e)+"_zy";

        TH2* proj = dynamic_cast<TH2*>(h->Project3D(projcmd.c_str()));
        proj->SetTitle(Form("#Delta #theta, CB, Element %d",e));
        proj->Draw("colz");

        proj->FitSlicesY(f,0,-1,0,"QNR");

        TH1D* h_sigma = nullptr;

        gDirectory->GetObject(Form("CB_theta_%d_zy_2",e), h_sigma);
        h_sigma->SetTitle(Form("#Delta #theta, fitted, CB, Element %d",e));

        TF1* thetaE = new TF1(Form("thetaE_%d",e),"expo+pol0(2)", 0, 1600);
        h_sigma->Fit(thetaE,"MRQ");

        TH1* p1 = proj->ProjectionY(Form("CB_theta_%d_proj",e),0,-1);
        TF1* p1f = new TF1("", "gaus", zbins.Start(), zbins.Stop());
        p1f->SetParameter(0, 1);
        p1f->SetParameter(1, 0);
        p1f->SetParameter(2,.1);
        p1->SetTitle(Form("#Delta #theta, global, CB, Element %d",e));
        p1->Fit(p1f,"QNR");

        cout << "Element " << e << ":  Global sigma="<< p1f->GetParameter(2) <<", E-Fitted: "
             << thetaE->GetParameter(0) << " "
             << thetaE->GetParameter(1) << " "
             << thetaE->GetParameter(2)
             << " Chi2/dof=" << thetaE->GetChisquare()/thetaE->GetNDF() << endl;

        h_global_sigma->SetBinContent(e+1,p1f->GetParameter(2));
        h_theta_fit_p0->SetBinContent(e+1, thetaE->GetParameter(0));
        h_theta_fit_p1->SetBinContent(e+1, thetaE->GetParameter(1));
        h_theta_fit_p2->SetBinContent(e+1, thetaE->GetParameter(2));

        fits->cd();
        if(first) {
            thetaE->Draw();
            first=false;
        } else {
            thetaE->Draw("same");
        }
        c->cd();


    }


    TH2CB* cb = new TH2CB("theta_sigmas", "#sigma #theta, global, CB");
    cb->SetZTitle("#sigma_{#theta} [rad]");
    cb->SetElements(*h_global_sigma);

    TH2CB* cba = new TH2CB("theta_sigma_p0", "#sigma #theta P0");
    cba->SetZTitle("p0");
    cba->SetElements(*h_theta_fit_p0);

    TH2CB* cbb = new TH2CB("theta_sigma_p1", "#sigma #theta P1");
    cbb->SetZTitle("p1");
    cbb->SetElements(*h_theta_fit_p1);

    TH2CB* cbc = new TH2CB("theta_sigma_p2", "#sigma #theta P2");
    cbc->SetZTitle("p2");
    cbc->SetElements(*h_theta_fit_p2);

    canvas("Sigma Theta CB")
            << cb
            << cba
            << cbb
            << cbc
            << endc;

}

void ExtractResolutions::AnalysePhiCB(TTree* tree)
{

    const unsigned nElements = 720;
    const BinSettings xbins(nElements, 0, nElements);
    const BinSettings ybins(16,     0.0, 1600.0);
    const BinSettings zbins(80,    -0.2,    0.2);


    TCanvas* c = new TCanvas();
    auto h = Draw(tree, "Phi:tE:Element", TCut(""), xbins, ybins, zbins);
    h->SetXTitle("Element");
    h->SetZTitle("#Delta #phi [rad]");
    h->SetYTitle("E [MeV]");
    h->SetName("CB_phi");

    TF1* f = new TF1("f", "gaus", zbins.Start(), zbins.Stop());
    f->SetParameter(0, 1);
    f->SetParameter(1, 0);
    f->SetParameter(2,.1);



    TCanvas* fits = new TCanvas();
    bool first = true;
    c->cd();


    TH1D* h_global_sigma = new TH1D("cb_phi_global_sigma", "Global #Delta #phi, cb", nElements, 0, nElements);
    TH1D* h_phi_fit_p0   = new TH1D("cb_phi_fit_p0", "Fit #sigma #phi, p0, cb", nElements, 0, nElements);
    TH1D* h_phi_fit_p1   = new TH1D("cb_phi_fit_p1", "Fit #sigma #phi, p1, cb", nElements, 0, nElements);
    TH1D* h_phi_fit_p2   = new TH1D("cb_phi_fit_p2", "Fit #sigma #phi, p2, cb", nElements, 0, nElements);

    for(int e=0; e<int(nElements); ++e) {
        h->GetXaxis()->SetRange(e,e+1);

        const string projcmd = to_string(e)+"_zy";

        TH2* proj = dynamic_cast<TH2*>(h->Project3D(projcmd.c_str()));
        proj->SetTitle(Form("#Delta #phi, CB, Element %d",e));
        proj->Draw("colz");

        proj->FitSlicesY(f,0,-1,0,"QNR");

        TH1D* h_sigma = nullptr;

        gDirectory->GetObject(Form("CB_phi_%d_zy_2",e), h_sigma);
        h_sigma->SetTitle(Form("#Delta #phi, fitted, CB, Element %d",e));

        TF1* phiE = new TF1(Form("phiE_%d",e),"expo+pol0(2)", 0, 1600);
        h_sigma->Fit(phiE,"MRQ");

        TH1* p1 = proj->ProjectionY(Form("CB_phi_%d_proj",e),0,-1);
        TF1* p1f = new TF1("", "gaus", zbins.Start(), zbins.Stop());
        p1f->SetParameter(0, 1);
        p1f->SetParameter(1, 0);
        p1f->SetParameter(2,.1);
        p1->SetTitle(Form("#Delta #phi, global, CB, Element %d",e));
        p1->Fit(p1f,"QNR");

        cout << "Element " << e << ":  Global sigma="<< p1f->GetParameter(2) <<", E-Fitted: "
             << phiE->GetParameter(0) << " "
             << phiE->GetParameter(1) << " "
             << phiE->GetParameter(2)
             << " Chi2/dof=" << phiE->GetChisquare()/phiE->GetNDF() << endl;

        h_global_sigma->SetBinContent(e+1,p1f->GetParameter(2));
        h_phi_fit_p0->SetBinContent(e+1, phiE->GetParameter(0));
        h_phi_fit_p1->SetBinContent(e+1, phiE->GetParameter(1));
        h_phi_fit_p2->SetBinContent(e+1, phiE->GetParameter(2));

        fits->cd();
        if(first) {
            phiE->Draw();
            first=false;
        } else {
            phiE->Draw("same");
        }
        c->cd();


    }


    TH2CB* cb = new TH2CB("phi_sigmas", "#sigma #phi, global, CB");
    cb->SetZTitle("#sigma_{#phi} [rad]");
    cb->SetElements(*h_global_sigma);

    TH2CB* cba = new TH2CB("phi_sigma_p0", "#sigma #phi P0");
    cba->SetZTitle("p0");
    cba->SetElements(*h_phi_fit_p0);

    TH2CB* cbb = new TH2CB("phi_sigma_p1", "#sigma #phi P1");
    cbb->SetZTitle("p1");
    cbb->SetElements(*h_phi_fit_p1);

    TH2CB* cbc = new TH2CB("phi_sigma_p2", "#sigma #phi P2");
    cbc->SetZTitle("p2");
    cbc->SetElements(*h_phi_fit_p2);

    canvas("Sigma Phi CB")
            << cb
            << cba
            << cbb
            << cbc
            << endc;

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
