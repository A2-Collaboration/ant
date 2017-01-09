#include "KinfitExtract.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TList.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "analysis/plot/root_draw.h"

#include "analysis/utils/Fitter.h"
#include "analysis/utils/uncertainties/Optimized.h"

#include <memory>

using namespace std;
using namespace ant;
using namespace ant::std_ext;


void KinfitExtract::Sweep1(TTree* t)
{
    const string p1name = "T_cb_photon_theta_const";
    const string p2name = "T_cb_photon_theta_Sin";
    const auto parm1vals = { 0.5, 2.0, 3.5 };
    const auto parm2vals = { 3.0, 5.5, 8.0 };

    auto g = new TGraph2D(int(parm1vals.size()*parm2vals.size()));
    const string gTitle = formatter() << "RMS d#theta;"<<p1name<<";"<<p2name<<";RMS";
    g->SetTitle(gTitle.c_str());

    canvas c;

    unsigned n = 0;
    for(const auto& p1 : parm1vals) {
        for(const auto& p2 : parm2vals) {

            const string title = formatter() << "#theta pull: const=" << p1 << ", Sin=" << p2;
            const string name  = formatter() << "_h_" << n;

            auto hist = new TH1D(name.c_str(), title.c_str(), 100, -3.0, 3.0);

            const string cmd   = formatter() << "1Pi0_Photon0_Theta_pull>>" << name;
            const string cut   = formatter() << "TaggW*(" << p1name << "==" << p1 << " && " << p2name << "=="  << p2 <<")";

            t->Draw(cmd.c_str(), cut.c_str());

            c << hist;

            const auto RMS = hist->GetRMS();

            g->SetPoint(int(n), p1, p2, RMS);

            ++n;
        }
    }

    c << endc;

    canvas("RMS") << drawoption("surf") << g << endc;

}

std::string KinfitExtract::JSONTest()
{
    analysis::utils::UncertaintyModels::Optimized_Oli1 m;
    return m.to_string_short();
}

string KinfitExtract::LoadAndDump(const string& s)
{
    analysis::utils::UncertaintyModels::Optimized_Oli1 m;
    m.load_from_string(s);
    return m.to_string();
}

void KinfitExtract::short_string_test()
{
    analysis::utils::UncertaintyModels::Optimized_Oli1 m;
    const auto s = m.to_string_simple();
    cout << s << endl;

    analysis::utils::UncertaintyModels::Optimized n;
    n.load_from_string_simple(s);
    cout << n.to_string_simple() << endl;
}

void AddMeanRMStoTitle(TH1* h, const std::string& ptitle) {
    const string title = formatter() << ptitle << "(Mean=" << h->GetMean() << ", RMS=" << h->GetRMS() << ")";
    h->SetTitle(title.c_str());
}

THStack* KinfitExtract::DrawPullSame(TH1* data, TH1* mc, const string& title, const bool draw)
{
//    std::pair<TH1*, TH1*> maxmin = {data, mc};

//    if(data->GetMaximum() < mc->GetMaximum())
//        swap(maxmin.first, maxmin.second);

 //   const double max = maxmin.first->GetMaximum();

    data->SetLineColor(kRed);
    data->Scale(1.0/data->Integral());
    mc->SetLineColor(kBlue);
    mc->Scale(1.0/mc->Integral());

    AddMeanRMStoTitle(data, "Data");
    AddMeanRMStoTitle(mc, "MC");

    auto h = new THStack();
    h->SetTitle(title.c_str());

    h->Add(data);
    h->Add(mc);

    if(draw) {
        auto c = new TCanvas();
        h->Draw("nostack");
        c->BuildLegend();
    }

    return h;

}

void KinfitExtract::DrawPullSameAll(TDirectory* data, TDirectory* mc)
{
    auto l = data->GetListOfKeys();

        if(!l)
            return;

        canvas stacks;

        TKey* key  = nullptr;
        TIter nextk(l);

        while((key = (TKey*)nextk()))
        {
            TH1* data_h = dynamic_cast<TH1*>(key->ReadObj());
            if ( !data_h )
                continue;

            TH1* mc_h = nullptr;
            mc->GetObject(data_h->GetName(), mc_h);

            if(mc_h) {
                auto stack = DrawPullSame(data_h, mc_h, data_h->GetTitle(), false);
                if(string_starts_with(string(data_h->GetName()), "p_")) {
                    stacks << drawoption("nostack") << padoption::Legend << stack;
                }
            }
        }

        stacks << endc;

}
