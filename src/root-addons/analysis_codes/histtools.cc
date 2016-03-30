#include "histtools.h"

#include <iostream>
#include <string>

#include "TH1.h"
#include "TH1D.h"

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
