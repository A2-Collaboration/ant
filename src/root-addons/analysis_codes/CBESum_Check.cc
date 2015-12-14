#include "CBESum_Check.h"

#include "TH2D.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "expconfig/ExpConfig.h"
#include "analysis/plot/root_draw.h"
#include "base/interval.h"
#include "base/std_ext/math.h"

#include <iostream>

using namespace std;
using namespace ant;

TH2CB* makeMeanHist(TH2D* h, interval<double> range) {
    h->GetXaxis()->SetRangeUser(range.Start(),range.Stop());

    TH2CB* cb = new TH2CB(std_ext::formatter() << "cb_" << h->GetName(),
                          std_ext::formatter() << "Mean " << range << " " << h->GetTitle() );
    ExpConfig::Setup::ManualName = "Setup_2014_07_EPT_Prod";
    auto cb_detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);


    interval<double> minmax(std_ext::NaN, std_ext::NaN);
    for(unsigned ch=0;ch<720;ch++) {
        if(cb_detector->IsIgnored(ch)) {
            cb->CreateMarker(ch);
            continue;
        }
        TH1D* h_proj = h->ProjectionX("h_proj",ch+1,ch+1);
        auto mean =  h_proj->GetMean();
        minmax.Extend(mean);
        cb->SetElement(ch,mean);
    }
    cb->GetZaxis()->SetRangeUser(minmax.Start(),minmax.Stop());
    return cb;
}

void CBESum_Check::AnalyseHist(TH2D* h_CBEsum, TH2D* h_E)
{
    auto cb_CBEsum = makeMeanHist(h_CBEsum, {400, 700});
    auto cb_E = makeMeanHist(h_E, {0, 1000});

    TH2CB* cb_div = dynamic_cast<TH2CB*>(cb_CBEsum->Clone());
    cb_div->Divide(cb_E);
    cb_div->GetZaxis()->UnZoom();

    canvas c;
    c << drawoption("colz")
      << padoption::enable(padoption::LogZ)
      << h_CBEsum
      << h_E
      << padoption::disable(padoption::LogZ)
      << cb_CBEsum
      << cb_E
      << cb_div
      << endc;
}

