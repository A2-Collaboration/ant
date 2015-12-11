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

void CBESum_Check::AnalyseHist(TH2D* h)
{
    h->GetXaxis()->SetRangeUser(400,700);

    TH2CB* cb = new TH2CB("cb","Mean CBEsum between 400 and 700 MeV");
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
    canvas c;
    c << drawoption("colz")
      << padoption::LogZ
      << h
      << drawoption("colz") << cb
      << samepad
      << drawoption("text") << new TH2CB()
      << endc;
}

