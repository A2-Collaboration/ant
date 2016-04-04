#include "ThetaCBToyMC.h"
#include "base/std_ext/math.h"
#include "TH2.h"
#include "TH3.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "histtools.h"
#include <iostream>
#include "TList.h"
#include "analysis/plot/root_draw.h"

using namespace std;
using namespace ant;
using namespace ant::MC;
using namespace ant::std_ext;

#include <cmath>
#include "TRandom2.h"

double ThetaCBToyMC::dTheta(const double theta, const double z, const double r) {
    const auto a0 = r*cos(theta);
    const auto b  = r*sin(theta);

    return atan2(b,a0-z);
}

TH2* ThetaCBToyMC::SimTheta(const unsigned n, const double r, const double l, TH2* hist)
{
    const double theta_min = degree_to_radian(20.0);
    const double theta_max = degree_to_radian(160.0);

    if(!hist) {
        hist = new TH2D("thetacb", Form("d#theta: %f cm target, %f cm radius", l, r), 180, radian_to_degree(theta_min), radian_to_degree(theta_max), 200,-20, 20);
        hist->SetXTitle("#theta [#circ]");
        hist->SetYTitle("d#theta [#circ]");
    }

    TRandom2 rng;
    rng.SetSeed();

    for(unsigned i=0; i<n; ++i) {
        const auto theta = rng.Uniform(theta_max - theta_min) + theta_min;
        const auto z     = rng.Uniform(l) - l/2.0;

        const auto dtheta = theta - dTheta(theta, z, r);

        hist->Fill(radian_to_degree(theta), radian_to_degree(dtheta));
    }

    return hist;

}

void ThetaCBToyMC::Analyse3(TH3* hist)
{
    const int EBins = hist->GetNbinsZ();

    const std::vector<int> colors = {kBlack, kRed, kBlue, kGreen, kPink, kCyan, kOrange, kTeal};

    new TCanvas();

    for(int iE=1; iE <= EBins; ++iE) {

        hist->GetZaxis()->SetRange(iE,iE);
        auto eSlice = dynamic_cast<TH2*>(hist->Project3D(Form("sliceE%d_yx", iE)));
        if(!eSlice) {
            cerr << "Projection failed!" << endl;
            return;
        }

        eSlice->SetTitle(Form("E > %f && E < %f", hist->GetZaxis()->GetBinLowEdge(iE), hist->GetZaxis()->GetBinUpEdge(iE)));

      //  ant::histtools::PlotMeansRMS(eSlice);

        auto g = ant::histtools::MeanRMS(eSlice);

        g.second->SetMarkerColor(colors.at( iE % colors.size()));
        g.second->SetMarkerStyle(7);

        if(iE==1)
            g.second->Draw("AP");
        else
            g.second->Draw("P same");


    }
}

void ThetaCBToyMC::SimThetaMulti(const unsigned n, const double rmin, const double rmax, const unsigned steps, const double l)
{
    const std::vector<int> colors = {kBlack, kRed, kBlue, kGreen, kPink, kCyan, kOrange, kTeal};

    canvas cHists("Thetas");
    auto cMean = new TCanvas("Mean");
    auto cRMS  = new TCanvas("RMS");

    for(unsigned i=0; i<steps; ++i) {
        const double r = rmin + ((rmax - rmin)/steps)*i;

        const auto hist = SimTheta(n, r, l);

        cHists << drawoption("colz")
                << hist;

        auto g = ant::histtools::MeanRMS(hist);

        g.first->SetMarkerStyle(7);
        g.second->SetMarkerStyle(7);
        g.first->SetMarkerColor(colors.at( i % colors.size()));
        g.second->SetMarkerColor(colors.at( i % colors.size()));

        g.first->SetTitle(Form("Mean r = %f", r));
        g.second->SetTitle(Form("RMS r = %f", r));

        if(i==0) {
            cMean->cd();
            g.first->Draw("AP");
            cRMS->cd();
            g.second->Draw("AP");
        } else {
            cMean->cd();
            g.first->Draw("P same");
            cRMS->cd();
            g.second->Draw("P same");
        }


    }

    cHists << endc;

}
