#include "TH1Ext.h"
#include "TH1.h"
#include "TH2.h"

using namespace ant;

void TH2Ext::MakeSameZRange(std::vector<TH2*> hists) {

    if(hists.empty())
        return;


    auto b1 = getZRange(*hists.front());

    for(auto it=next(hists.begin()); it!=hists.end(); ++it) {
        b1.Extend(getZRange(**it));
    }

    for(auto h : hists) {
        h->SetMinimum(b1.Start());
        h->SetMaximum(b1.Stop());
    }

}

void TH2Ext::ClearHistogram(TH2D* hist, const double v) {
    for(int x=1; x<=hist->GetNbinsX(); ++x) {
        for(int y=1; y<=hist->GetNbinsY(); ++y) {
            hist->SetBinContent(x,y, v);
        }
    }
}

interval<double> TH2Ext::getZRange(const TH2& hist) {
    return {hist.GetMinimum(), hist.GetMaximum()};
}

double ant::GetMaxPos(TH1* hist){
    return hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
}
