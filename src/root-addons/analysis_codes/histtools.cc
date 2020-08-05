#include "histtools.h"

#include "base/std_ext/string.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TVirtualX.h"

using namespace std;
using namespace ant;

TH1* histtools::CutPosPlot(const TH1* const h1, const TH1* const h2)
{
    if(h1->GetNbinsX() != h2->GetNbinsX()) {
        cerr << "ERROR: Bin count mismatch!" << endl;
    }

    auto h = new TH1D("","",h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
    h->SetXTitle("Cut Pos");
    h->SetYTitle(string("#frac{" + string(h1->GetTitle()) + "}{" + string(h2->GetTitle()) + "}").c_str());

    if(!h)
        cerr << "Error cloning histogram!" << endl;

    double sum1 = 0.0;
    double sum2 = 0.0;
    for(int i=1; i<=h1->GetNbinsX(); ++i) {
        sum1 += h1->GetBinContent(i);
        sum2 += h2->GetBinContent(i);

        const auto v = sum1 / sum2;

        if(isfinite(v))
            h->SetBinContent(i, v);

    }
    return h;
}

TH1*histtools::PlotSignificance(const TH1* const h1, const TH1* const h2)
{
    if(h1->GetNbinsX() != h2->GetNbinsX()) {
        cerr << "ERROR: Bin count mismatch!" << endl;
    }

    auto h = new TH1D("","",h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
    h->SetXTitle("Cut Pos");
    h->SetYTitle("Significance");

    if(!h)
        cerr << "Error cloning histogram!" << endl;

    double sum1 = 0.0;
    double sum2 = 0.0;
    for(int i=1; i<=h1->GetNbinsX(); ++i) {
        sum1 += h1->GetBinContent(i);
        sum2 += h2->GetBinContent(i);

        const auto v = sum1 / sqrt(sum1+sum2);

        if(isfinite(v))
            h->SetBinContent(i, v);

    }
    return h;
}

std::pair<TGraph*, TGraph*> histtools::MeanRMS(const TH2* const h)
{
    auto gmeans = new TGraph(h->GetNbinsX());
    gmeans->SetTitle(Form("Y-Means: %s", h->GetTitle()));
    gmeans->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    gmeans->GetYaxis()->SetTitle(Form("Mean(%s)", h->GetYaxis()->GetTitle()));
    gmeans->SetMarkerStyle(6);


    auto grms   = new TGraph(h->GetNbinsX());
    grms->SetTitle(Form("Y-RMS: %s", h->GetTitle()));
    grms->SetMarkerStyle(6);
    grms->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    grms->GetYaxis()->SetTitle(Form("RMS(%s)", h->GetYaxis()->GetTitle()));

    for(int i=1; i<h->GetNbinsX(); ++i) {
        auto proj = h->ProjectionY("",i,i);
        const auto mean = proj->GetMean();
        const auto rms  = proj->GetRMS();
        const auto x    = h->GetXaxis()->GetBinCenter(i);

        gmeans->SetPoint(i, x, mean);
        grms->SetPoint(i, x, rms);
        delete proj;
    }

    return {gmeans, grms};
}

void histtools::PlotMeansRMS(const TH2* const h)
{
    auto graphs = MeanRMS(h);

    auto c = new TCanvas();
    c->SetTitle(h->GetTitle());

    c->Divide(2,1);
    c->cd(1);
    graphs.first->Draw("AP");
    c->cd(2);
    graphs.second->Draw("AP");

}

struct avg_t {
    double sum = 0.0;
    long long n = 0;
    avg_t& operator()(double v) {
        if(isfinite(v))
        { sum+=v; ++n; }
        return *this;
    }

    double GetAverage() const { return sum / n; }
    double GetSum() const { return sum; }
};

TH1D *histtools::ProjectX(TH2 *hist)
{
    const auto xbins = hist->GetNbinsX();
    const auto ybins = hist->GetNbinsY();
    TH1D* res = new TH1D("","", xbins, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    for(int x = 1; x<= xbins; ++x) {
        avg_t f;
        for(int y = 1; y<= ybins; ++y) {
            f(hist->GetBinContent(x,y));
        }

        res->SetBinContent(x, hist->TestBit(TH1::kIsAverage) ? f.GetAverage() : f.GetSum());
    }

    return res;
}

void histtools::CleanupLabeledHist(TH1* h)
{
    if(h->GetXaxis()->GetLabels() == nullptr) {
        cerr << "Can only work on labeled histograms" << endl;
        return;
    }

    auto h_clean = new TH1D((string(h->GetName())+"_clean").c_str(),h->GetTitle(), 1, 0, 1);
    h_clean->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
    h_clean->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());

    vector<pair<string, double>> contents;
    for(int bin=1;bin<h->GetNbinsX();bin++) {
        const auto entries =  h->GetBinContent(bin);
        if(entries<1)
            continue;
        contents.emplace_back(make_pair(h->GetXaxis()->GetBinLabel(bin), entries));
    }

    sort(contents.begin(), contents.end(), [] (const pair<string, double>& a, const pair<string, double>& b) {
        return a.second > b.second;
    });

    for(auto& i : contents) {
        h_clean->Fill(i.first.c_str(), i.second);
    }

    new TCanvas;
    h_clean->Draw();
    h_clean->GetXaxis()->SetRangeUser(0,contents.size());
    gPad->SetLogy();
    {
        int bin = 1;
        for(auto& i : contents) {
            const auto entries = i.second;
            auto lbl = new TLatex(bin-0.5, entries*1.1, string(std_ext::formatter() << entries).c_str());
            lbl->SetTextAngle(90);
            lbl->Draw();
            bin++;
        }
    }
}
