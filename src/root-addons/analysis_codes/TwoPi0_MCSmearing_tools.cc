#include "TwoPi0_MCSmearing_tools.h"

#include "base/TH_ext.h"
#include "root-addons/analysis_codes/Math.h"
#include "base/std_ext/math.h"
#include "root-addons/cbtaps_display/TH2CB.h"
#include "analysis/plot/root_draw.h"
#include "base/std_ext/string.h"
#include "base/BinSettings.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TList.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TAxis.h"
#include "TROOT.h"

#include <list>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;


struct HAxis_t {
    std::string Title;
    BinSettings bins;
    HAxis_t(const TAxis* axis): Title(axis->GetTitle()), bins(unsigned(axis->GetNbins()), axis->GetXmin(), axis->GetXmax()) {}
};


class InspectorCanvas : public TCanvas {
protected:
    int obx = -1;
    int oby = -1;
    const string histogram_base_name;
    std::list< pair<string, TDirectory*> > dirs;

public:
    InspectorCanvas(const string& TH1_base_name, list<pair<string, TDirectory*>> directories): TCanvas(), histogram_base_name(TH1_base_name), dirs(directories) {}

    virtual ~InspectorCanvas() {}

    virtual void HandleInput(const EEventType button, const Int_t px, const Int_t py) override;
};

void InspectorCanvas::HandleInput(const EEventType button, const Int_t px, const Int_t py)
{
    if( button == kButton2Motion || button == kButton1Up) {

        auto h2 = dynamic_cast<TH2*>(fClickSelected);

        if(h2) {

            Double_t x,y;
            AbsPixeltoXY(px,py,x,y);

            const auto bin = h2->FindBin(x,y);
            int bx ,by, dummy;
            h2->GetBinXYZ(bin, bx, by, dummy);

            if(bx != obx || by != oby) {

                const string hist_name = formatter() << histogram_base_name << "_" << bx << "_" << by;

                for(auto d : dirs) {
                    TH1* h = nullptr;
                    d.second->GetObject(hist_name.c_str(), h);

                    if(h) {

                        const string cname = formatter() << size_t(this) << "_" << d.first;
                        auto obj = gROOT->FindObjectAny(cname.c_str());
                        TCanvas* c = dynamic_cast<TCanvas*>(obj);
                        if(!c) {
                            c = new TCanvas(cname.c_str(), d.first.c_str());
                        }

                        c->cd();
                        gStyle->SetOptFit(1);
                        h->Draw();
                        c->Modified();
                        c->Update();
                        obx = bx;
                        oby = by;

                    } else {
                        cout << "TH1 " << hist_name << " not found" << endl;
                    }
                }
            }
        }
    }

    TCanvas::HandleInput(button, px, py);
}

TCanvas* TwoPi0_MCSmearing_Tool::getInspectorCanvas(TH2* h, const string& hist_base, TDirectory* dir, const string& n1)
{
    if(!dir)
        dir = gDirectory;

    auto c = new InspectorCanvas(hist_base, {{n1, dir}});

    h->Draw("colz");
    c->SetEditable(false);

    return c;
}

TCanvas* TwoPi0_MCSmearing_Tool::getInspectorCanvas(TH2* h, const string& hist_base, TDirectory* dir1, const string& n1, TDirectory* dir2, const string& n2)
{

    auto c = new InspectorCanvas(hist_base, {{n1, dir1}, {n2, dir2}});

    h->Draw("colz");
    c->SetEditable(false);

    return c;
}

PeakFitResult_t ant::TwoPi0_MCSmearing_Tool::Fit(TH1* h, const std::string& prefix, const bool verbose)
{
    if(h->GetEntries() <= 10000)
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

void TwoPi0_MCSmearing_Tool::DrawSame(TH1* h1, TH1* h2)
{
    const auto maxy = max(h1->GetMaximum(), h2->GetMaximum());

    h1->GetYaxis()->SetRangeUser(0,maxy);
    h2->GetYaxis()->SetRangeUser(0,maxy);

    h1->Draw();
    h2->Draw("same");

    gPad->BuildLegend();

}

channelFitResult_t TwoPi0_MCSmearing_Tool::FitChannel(TH2* h2, const int ch)
{
    auto h = h2->ProjectionX(Form("proj_%d",ch), ch, ch+1);
    return {Fit(h,h->GetName()), h};
}

MultiChannelFitResult_t TwoPi0_MCSmearing_Tool::FitAllChannels(TH2* h2)
{
    const auto ElementAxis = HAxis_t(h2->GetYaxis());

    auto g_sigmas = new TH1D("sigma", Form("%s #simga", h2->GetName()), ElementAxis.bins.Bins(),ElementAxis.bins.Start(), ElementAxis.bins.Stop() );
    g_sigmas->GetXaxis()->SetTitle(ElementAxis.Title.c_str());
    g_sigmas->GetYaxis()->SetTitle("#sigma [MeV]");

    auto g_pos    = new TH1D("pos",   Form("%s pos",    h2->GetName()), ElementAxis.bins.Bins(),ElementAxis.bins.Start(), ElementAxis.bins.Stop() );
    g_pos->GetXaxis()->SetTitle(ElementAxis.Title.c_str());
    g_pos->GetYaxis()->SetTitle("Peak Pos [MeV]");

    for(int i=1; i <= int(ElementAxis.bins.Bins()); ++i) {
        cout << i << "\r" << flush;
        auto r = FitChannel(h2, i);
        const auto value_ok = isfinite(r.result.chi2dof);
        if(value_ok) {
            g_sigmas->SetBinContent(i, r.result.sigma );
            g_pos->SetBinContent   (i,  r.result.pos  );
        }

    }

    auto c =new TCanvas();
    c->Divide(2,2);
    c->cd(1);
    g_sigmas->Draw("P");
    c->cd(2);
    c->cd(3);
    g_pos->Draw("P");


    return {g_pos, g_sigmas};
}

TH1* THDataMCDiff(const TH1* mc, const TH1* data, const string& name) {

    auto diff = TH_ext::Clone(data, name);
    diff->Reset();

    for(int i=0;i<=data->GetNbinsX(); ++i) {

        const auto d = data->GetBinContent(i);
        const auto m = mc->GetBinContent(i);

        if(isfinite(d) && isfinite(m)) {
            diff->SetBinContent(i, (d-m)/m);
        }
    }

    return diff;
}

TH2* THDataMCDiff(const TH2* mc, const TH2* data, const string& name) {

    auto diff = TH_ext::Clone(data, name);
    diff->Reset();

    double min= inf;
    double max=-inf;
    for(int x=0;x<=data->GetNbinsX(); ++x) {
        for(int y=0;y<=data->GetNbinsY(); ++y) {

            const auto d = data->GetBinContent(x,y);
            const auto m = mc->GetBinContent(x,y);

            if(isfinite(d) && d > 0.0 && isfinite(m) && m > 0.0) {
                const auto v = (d-m)/m;
                min=std::min(min,v);
                max=std::max(max,v);
                diff->SetBinContent(x, y, v);
            } else {
                diff->SetBinContent(x,y, 0.0);
            }
        }
    }
    diff->SetTitle("(data-mc)/mc");

    //diff->GetZaxis()->SetRangeUser(min,max);

    return diff;
}

//TODO: Make work for EElement too!
void TwoPi0_MCSmearing_Tool::CompareMCData2D(TDirectory* mc, TDirectory* data, const string& folder)
{
    TH2* mc_sigma = nullptr;
    mc->GetObject((folder+"/sigma").c_str(), mc_sigma);
    if(!mc_sigma) {
        cerr << "MC hist not found!" << endl;
        return;
    }

    TH2* mc_pos = nullptr;
    mc->GetObject((folder+"/pos").c_str(), mc_pos);
    if(!mc_pos) {
        cerr << "MC hist not found!" << endl;
        return;
    }

    TH2* data_sigma = nullptr;
    data->GetObject((folder+"/sigma").c_str(), data_sigma);
    if(!data_sigma) {
        cerr << "data hist not found!" << endl;
        return;
    }

    TH2* data_pos = nullptr;
    data->GetObject((folder+"/pos").c_str(), data_pos);
    if(!data_pos) {
        cerr << "data hist not found!" << endl;
        return;
    }

    auto datamc_sigma = THDataMCDiff(mc_sigma, data_sigma, "datamc_simga");
    datamc_sigma->GetZaxis()->SetRangeUser(-.2,.6);

    auto datamc_pos = THDataMCDiff(mc_pos, data_pos, "datamc_pos");
    datamc_pos->GetZaxis()->SetRangeUser(-1111,-1111);

    auto sigma_s = TH_ext::Apply(mc_sigma, data_sigma,
                                 [] (const double& mc, const double& data) -> double {
        const auto v = sqrt(sqr(data) - sqr(mc));
        return isfinite(v) ? v : 0.0;

    });

    sigma_s->SetName("energy_initial");
    sigma_s->SetTitle("Energy Smearing [MeV] (only use for initial step)");

    canvas("Data/MC") << drawoption("colz") << datamc_sigma << datamc_pos << sigma_s << endc;

}

TH2* TwoPi0_MCSmearing_Tool::CalculateInitialSmearing(const TH2* sigma_data, const TH2* sigma_MC)
{

    auto d = TH_ext::Apply(sigma_data, sigma_MC,
                                 [] (const double& data, const double& mc) -> double {
        const auto v = data - mc;
        return isfinite(v) ? v : -1.0;

    });
    d->SetName("im_d");

    const auto maxdiff = d->GetMaximum();

    auto sigma_s = TH_ext::Apply(sigma_data, sigma_MC,
                                 [maxdiff] (const double& , const double& ) -> double {
        return maxdiff;

    });
    sigma_s->SetTitle("Energy smearing");
    sigma_s->SetName("energy_smearing");



    return sigma_s;
}

TH2* TwoPi0_MCSmearing_Tool::CalculateUpdatedSmearing(const TH2* sigma_data, const TH2* current_sigma_MC, const TH2* last_smear)
{

    const auto d =  TH_ext::Apply(sigma_data, current_sigma_MC, [] (const double d, const double mc) { return d-mc;});
    d->SetName("im_d");

    const auto hs = vector<const TH2*>({sigma_data, current_sigma_MC, last_smear});
    auto factor = TH_ext::ApplyMany(hs, [] (const vector<double>& v) {
        const auto s = v.at(2) + v.at(2) * 0.5 * (v.at(0)/v.at(1)-1);
        return max(s, 0.0);
    });
    factor->SetName("energy_smearing");
    factor->SetTitle("Energy smearing");

    canvas("Step") << drawoption("colz") << factor << endc;

    return factor;

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


void TwoPi0_MCSmearing_Tool::CompareMCData1D(TDirectory* mc, TDirectory* data)
{
    TH1* mc_sigma = nullptr;
    mc->GetObject("Element/cb_pi0_h_sigma", mc_sigma);
    if(!mc_sigma) {
        cerr << "MC hist not found!" << endl;
        return;
    }

    TH1* mc_pos = nullptr;
    mc->GetObject("Element/cb_pi0_h_pos", mc_pos);
    if(!mc_pos) {
        cerr << "MC hist not found!" << endl;
        return;
    }

    TH1* data_sigma = nullptr;
    data->GetObject("Element/cb_pi0_h_sigma", data_sigma);
    if(!data_sigma) {
        cerr << "data hist not found!" << endl;
        return;
    }

    TH1* data_pos = nullptr;
    data->GetObject("Element/cb_pi0_h_pos", data_pos);
    if(!data_pos) {
        cerr << "data hist not found!" << endl;
        return;
    }

    auto mg = new THStack();
    mc_sigma->SetMarkerColor(kRed);
    mc_sigma->SetMarkerStyle(7);
    mc_sigma->SetTitle(Form("MC: %s", mc_sigma->GetTitle()));
    mg->Add(mc_sigma);
    data_sigma->SetMarkerColor(kBlue);
    data_sigma->SetMarkerStyle(7);
    data_sigma->SetTitle(Form("Data: %s", data_sigma->GetTitle()));
    mg->Add(data_sigma);


    auto datamcdiff_sigma = THDataMCDiff(mc_sigma, data_sigma, "datamc_sigma");
    datamcdiff_sigma->SetTitle("#sigma: (data-mc)/data");
    datamcdiff_sigma->SetMarkerStyle(7);

    auto datamcdiff_pos = THDataMCDiff(mc_pos, data_pos, "datamc_pos");
    datamcdiff_pos->SetTitle("#sigma: (data-mc)/data");
    datamcdiff_pos->SetMarkerStyle(7);

    auto c = new TCanvas();
    c->Divide(2,2);

    c->cd(1);
    mg->Draw("nostack P");
    mg->GetXaxis()->SetTitle(mc_sigma->GetXaxis()->GetTitle());
    mg->GetYaxis()->SetTitle(mc_sigma->GetYaxis()->GetTitle());
    gPad->BuildLegend();

    c->cd(2);
    datamcdiff_sigma->Draw("P");

    c->cd(3);
    datamcdiff_pos->Draw("P");

    c->cd(4);

    auto h_cb = new TH2CB("cb_data_mcdiff","(Data-mc)/data #sigma");
    h_cb->SetElements(*datamcdiff_sigma);
    h_cb->Draw("colz");
}


TH2D* makeHist(const string& name, const string& title, const HAxis_t& xaxis, const HAxis_t& yaxis) {
    auto h = new TH2D(name.c_str(), title.c_str(), int(xaxis.bins.Bins()), xaxis.bins.Start(), xaxis.bins.Stop(), int(yaxis.bins.Bins()), yaxis.bins.Start(), yaxis.bins.Stop());
    h->SetXTitle(xaxis.Title.c_str());
    h->SetYTitle(yaxis.Title.c_str());
    return h;
}

TH2*TwoPi0_MCSmearing_Tool::AnalyseChannelE(TH3* h3)
{
    const HAxis_t channelBins(h3->GetZaxis());
    const HAxis_t eBins(h3->GetYaxis());

    TH2* h_sigmas = makeHist("sigma", "Sigma",            eBins, channelBins);
    TH2* h_pos    = makeHist("pos", "#pi^{0} peak pos", eBins, channelBins);
    TH2* stat     = makeHist("stats", "Statistics",       eBins, channelBins);
    TH2* status   = makeHist("status", "Fit status",       eBins, channelBins);
    TH2* chi2dof  = makeHist("chi2dof", "Chi2/dof",         eBins, channelBins);

    for(int element=1;element<=int(channelBins.bins.Bins());++element) {
        for(int ebin=1;ebin<=int(eBins.bins.Bins());++ebin) {
            cout << "-> " << setfill('0') << setw(2) << element << ":" <<  setfill('0') << setw(2) << ebin << "\r" << flush;
            auto h = h3->ProjectionX(Form("proj_%d_%d",ebin,element),ebin,ebin+1,element,element+1,"");
            const auto r = Fit(h,h->GetName());

            stat->SetBinContent(ebin, element, h->GetEntries());
            status->SetBinContent(ebin, element, r.status);

            const auto value_ok = (r.status==4000 && isfinite(r.chi2dof));
            if(value_ok) {
                h_sigmas->SetBinContent(ebin, element, abs(r.sigma) );
                h_pos->SetBinContent   (ebin, element, r.pos );
                chi2dof->SetBinContent (ebin, element, r.chi2dof );
            }

        }
    }

    h_pos->GetZaxis()->SetRangeUser(130,140);

    canvas c("Per Element Fits");
    c << drawoption("colz") << h_sigmas << chi2dof << stat << h_pos << endc;

    return h_sigmas;
}

TH2* TwoPi0_MCSmearing_Tool::RelativeDiff(const TH2* h1, const TH2* h2)
{
    auto res = dynamic_cast<TH2*>(h1->Clone());

    res->Add(h2,-1.0);

    res->Divide(h2);

    return res;
}






TBinGraphCanvas::TBinGraphCanvas(): TCanvas(), graph(nullptr), obx(-1), oby(-1)
{

}

TBinGraphCanvas::~TBinGraphCanvas()
{
    delete graph;
}

void TBinGraphCanvas::Add(TH2* h)
{
    if(hists.empty() || TH_ext::haveSameBinning(hists.front(), h)) {
        hists.push_back(h);
        this->cd();
        h->Draw("colz");
        if(graph) {
            delete graph;
        }

        graph = new TGraph(int(hists.size()));
    }
}

void TBinGraphCanvas::HandleInput(const EEventType button, const Int_t px, const Int_t py)
{
    if( button == kButton2Motion || button == kButton1Up) {

        auto h2 = dynamic_cast<TH2*>(fClickSelected);

        if(h2) {

            Double_t x,y;
            AbsPixeltoXY(px,py,x,y);

            const auto bin = h2->FindBin(x,y);
            int bx ,by, dummy;
            h2->GetBinXYZ(bin, bx, by, dummy);

            if(bx != obx || by != oby) {

                FillGraph(bx, by);


                const string cname = formatter() << size_t(this);
                auto obj = gROOT->FindObjectAny(cname.c_str());
                TCanvas* c = dynamic_cast<TCanvas*>(obj);
                if(!c) {
                    c = new TCanvas(cname.c_str(), "");
                }

                c->cd();
                gStyle->SetOptFit(1);
                graph->Draw("APL");
                c->Modified();
                c->Update();
                obx = bx;
                oby = by;


            }
        }
    }
    TCanvas::HandleInput(button, px, py);
}

void TBinGraphCanvas::FillGraph(const int x, const int y)
{
    int i=0;
    for(const auto& h : hists) {
        graph->SetPoint(i, double(i), h->GetBinContent(x,y));
        ++i;
    }
}

