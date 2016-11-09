#include "TwoPi0_MCSmearing_tools.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TList.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include <cmath>
#include "TGraph.h"
#include "TDirectory.h"
#include "root-addons/analysis_codes/Math.h"
#include "base/std_ext/math.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TAxis.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include <iostream>
#include <iomanip>
#include "analysis/plot/root_draw.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;


PeakFitResult_t ant::TowPi0_MCSmearing_Tool::Fit(TH1* h, const std::string& prefix, const bool verbose)
{
    if(h->GetEntries() <= 0)
        return PeakFitResult_t(NaN, NaN, NaN, -10);

    const double r_min = 80.0;
    const double r_max = 200.0;
    const double exp_center = 135.0;
    const double exp_width= 12.0;
    const int npx= 300;

    TF1* sig = new TF1((prefix+string("_sig")).c_str(), "gaus", r_min, r_max);
    sig->SetLineColor(kGreen);
    sig->SetNpx(npx);

    // height
    sig->SetParameter(0, 0.5 * h->GetMaximum());

    sig->SetParameter(1, exp_center);

    // width
    sig->SetParameter(2, 12.0);



    TF1* bg = new TF1((prefix+string("_bg")).c_str(), "pol3", r_min, r_max);
    bg->SetLineColor(kBlue);

    bg->SetParameter(0,0);
    bg->SetParName(0, "BG p_{0}");
    bg->SetParameter(1,0);
    bg->SetParName(1, "BG p_{1}");
    bg->SetParameter(2,0);
    bg->SetParName(2, "BG p_{2}");
    bg->SetParameter(3,0);
    bg->SetParName(3, "BG p_{3}");

    TFSum::FitRanged(h, bg, r_min, exp_center - 2*exp_width, exp_center+2*exp_width, r_max, "REM0NBQ");
    //bg->FixParameter(0, bg->GetParameter(0));
    //bg->FixParameter(1, bg->GetParameter(1));
    //bg->FixParameter(2, bg->GetParameter(2));


    TFSum* sum = new TFSum((prefix+string("_sum")).c_str(), sig, bg, r_min, r_max);
    sum->SetNpx(npx);

    h->SetStats(true);
    if(verbose) {
        TCanvas* c = new TCanvas();
        c->SetTitle(Form("Fit to %s", h->GetTitle()));


        gStyle->SetOptFit(1);
        h->Draw();
    }

    auto fitopt = [] (const bool& verbose) -> string { if(verbose) return "SREM0NB"; else return "SREM0NBQ"; };

    const auto fr = h->Fit(sum->Function(), fitopt(verbose).c_str());


    sum->SyncToFcts();

    if(verbose)
        sum->Draw();

    h->GetListOfFunctions()->Add(sum->Function());
    h->GetListOfFunctions()->Add(bg);

    const auto peak_pos  = sig->GetParameter(1);

    return PeakFitResult_t(fr->Chi2()/fr->Ndf(), peak_pos, sig->GetParameter(2), fr->Status());

}

void TowPi0_MCSmearing_Tool::DrawSame(TH1* h1, TH1* h2)
{
    const auto maxy = max(h1->GetMaximum(), h2->GetMaximum());

    h1->GetYaxis()->SetRangeUser(0,maxy);
    h2->GetYaxis()->SetRangeUser(0,maxy);

    h1->Draw();
    h2->Draw("same");

    gPad->BuildLegend();

}

channelFitResult_t TowPi0_MCSmearing_Tool::FitChannel(TH2* h2, const int ch)
{
    auto h = h2->ProjectionX(Form("%s_ch%d",h2->GetName(),ch), ch, ch+1);
    return {Fit(h,h->GetName()), h};
}

MultiChannelFitResult_t TowPi0_MCSmearing_Tool::FitAllChannels(TH2* h2, const std::string& prefix)
{
    TGraph* g = new TGraph(h2->GetNbinsY());
    g->SetName((prefix+string("_sigmas_g")).c_str());
    g->SetTitle((prefix+string(" #sigma")).c_str());
   // gDirectory->Add(g);
    TH1D*   h = new TH1D((prefix+string("_sigmas_h")).c_str(),(prefix+string(" #sigma")).c_str(), 30, 0, 30);

    int p=0;
    for(int i=0; i<h2->GetNbinsY(); ++i) {
        auto r = FitChannel(h2, i);
        if(isfinite(r.result.chi2dof)) {
            g->SetPoint(p++, i, r.result.sigma);
            h->Fill(r.result.sigma);
        }
    }

    auto c =new TCanvas();
    c->Divide(2,1);
    c->cd(1);
    g->Draw("AP");
    g->GetXaxis()->SetTitle("Element Nr.");
    g->GetYaxis()->SetTitle("#sigma [MeV]");
    c->cd(2);
    h->SetXTitle("#sigma [MeV]");
    h->Draw();

    return {g, h};
}

inline TGraph* DiffGraph(TGraph* g1, TGraph* g2) {
    if(g1->GetN() != g2->GetN()) {
        cerr << "GetN() does not match" << endl;
        return nullptr;
    }

    TGraph* r = new TGraph(g1->GetN());
    for(int i=0;i<g1->GetN(); ++i) {
        double x1, y1;
        g1->GetPoint(i,x1,y1);
        double x2, y2;
        g2->GetPoint(i,x2,y2);

        if(x1 != x2) {
            cerr << "Point " << i << "x does not match" << endl;
            delete r;
            return nullptr;
        }

        r->SetPoint(i,x1, (y1-y2)/y1);
    }

    return r;
}

inline void ProjectGraph(TH1* h, const TGraph* g) {
    for(int i=0;i<g->GetN(); ++i) {
        double dummy, y;
        g->GetPoint(i,dummy,y);
        h->Fill(y);
    }
}

inline TH2CB* ProjectGraphCB(const TGraph* g) {
    auto cb = new TH2CB(Form("h_%s",g->GetName()), "");

    if(g->GetN() != cb->GetNumberOfElements()) {
        cerr << "Grapth does not have " << cb->GetNumberOfElements() << " entries!" << endl;
        return cb;
    }

    for(int i=0; i<cb->GetNumberOfElements(); ++i) {
        double x,y;
        g->GetPoint(i,x,y);
        cb->SetElement(unsigned(x),y);
    }

    return cb;
}

void TowPi0_MCSmearing_Tool::CompareMCData(TDirectory* mc, TDirectory* data)
{
    TH2* mc_h2 = nullptr;
    mc->GetObject("TwoPi0_MCSmearing/m2Pi0/cb_pi0", mc_h2);
    if(!mc_h2) {
        cerr << "MC hist not found!" << endl;
        return;
    }

    TH2* data_h2 = nullptr;
    data->GetObject("TwoPi0_MCSmearing/m2Pi0/cb_pi0", data_h2);
    if(!data_h2) {
        cerr << "Data hist not found!" << endl;
        return;
    }

    auto dir = gDirectory;
    auto mc_Dir = dir->mkdir("mc_hists","");
    auto data_Dir = dir->mkdir("data_hists");
    mc_Dir->cd();
    auto mc_res = FitAllChannels(mc_h2, "mc");
    data_Dir->cd();
    auto data_res = FitAllChannels(data_h2, "data");
    dir->cd();

    TMultiGraph* mg = new TMultiGraph();
    mc_res.g->SetMarkerColor(kRed);
    mc_res.g->SetMarkerStyle(7);
    mg->Add(mc_res.g);
    data_res.g->SetMarkerColor(kBlue);
    data_res.g->SetMarkerStyle(7);
    mg->Add(data_res.g);

    THStack* s = new THStack();
    mc_res.h->SetLineColor(kRed);
    s->Add(mc_res.h);
    data_res.h->SetLineColor(kBlue);
    s->Add(data_res.h);

    auto datamcdiff = DiffGraph(data_res.g, mc_res.g);
    datamcdiff->SetMarkerStyle(7);

    TH1D* h_datamcdiff = new TH1D("h_datamcdiff", "", 40, -.4, .4);
    ProjectGraph(h_datamcdiff, datamcdiff);

    auto c = new TCanvas();
    c->Divide(2,2);
    c->cd(1);
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle(mc_res.g->GetXaxis()->GetTitle());
    mg->GetYaxis()->SetTitle(mc_res.g->GetYaxis()->GetTitle());
    gPad->BuildLegend();
    c->cd(3);
    datamcdiff->Draw("AP");
    datamcdiff->GetXaxis()->SetTitle(mc_res.g->GetXaxis()->GetTitle());
    datamcdiff->GetXaxis()->SetTitle("(data - mc), relative");
    c->cd(2);
    s->Draw("nostack");
    s->GetXaxis()->SetTitle(mc_res.h->GetXaxis()->GetTitle());
    s->GetYaxis()->SetTitle(mc_res.h->GetYaxis()->GetTitle());
    gPad->BuildLegend();
    c->cd(4);
    h_datamcdiff->Draw();

    auto cb_canvas = new TCanvas();
    auto h_cb = ProjectGraphCB(datamcdiff);
    h_cb->SetTitle("(Data-mc)/data #sigma");
    h_cb->Draw("colz");
}

TH2*TowPi0_MCSmearing_Tool::AnalyseChannelE(TH3* h3)
{
    const auto channels = h3->GetNbinsZ();
    const auto ebins    = h3->GetNbinsY();
    const auto Emax     = h3->GetYaxis()->GetXmax();
    TH2* res    = new TH2D(Form("ch_e_%s",       h3->GetName()),"Sigma",ebins,0,Emax,channels,0,channels );
    TH2* stat   = new TH2D(Form("ch_e_%s_stats", h3->GetName()),"Statistics",ebins,0,Emax,channels,0,channels );
    TH2* status = new TH2D(Form("ch_e_%s_status", h3->GetName()),"Fit statys",ebins,0,Emax,channels,0,channels );
    TH2* chi2dof = new TH2D(Form("ch_e_%s_chi2dof", h3->GetName()),"Chi2/dof",ebins,0,Emax,channels,0,channels );

    for(int element=1;element<=70;++element) {
        for(int ebin=1;ebin<=ebins;++ebin) {
            cout << "-> " << setfill('0') << setw(2) << element << ":" <<  setfill('0') << setw(2) << ebin << "\r" << flush;
            auto h = h3->ProjectionX(Form("p_%s_%d_%d",h3->GetName(),ebin,element),ebin,ebin+1,element,element+1,"");
            const auto r = Fit(h,h->GetName());

            stat->SetBinContent(ebin, element, h->GetEntries());
            status->SetBinContent(ebin, element, r.status);

            if(isfinite(r.status==4000)) {
                res->SetBinContent(ebin, element, r.sigma);
                chi2dof->SetBinContent(ebin, element, r.chi2dof);
            }
        }
    }

    canvas c("Per Element Fits");
    c << drawoption("colz") << res << chi2dof << stat << status << endc;

    return res;
}
