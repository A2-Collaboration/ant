#include "OmegaEtaG.h"

#include "TTree.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TCut.h"
#include "THStack.h"

#include "analysis/plot/root_draw.h"
#include "analysis/plot/Histogram.h"

#include "TreeTools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;

void OmegaEtaG::Chi2Plot(TTree* tree)
{
    new TCanvas("chi2_plot","Kinfit Chi2");
    TH1* chi2_bg = Draw(tree, "EPB_chi2dof","SigBgFlag==2", 1000,0,10);
    chi2_bg->SetTitle("#chi^{2} Background");
    chi2_bg->SetLineColor(kRed);

    TH1* chi2_sig = Draw(tree, "EPB_chi2dof","SigBgFlag!=2", 1000,0,10);
    chi2_sig->SetTitle("#chi^{2} Signal + Reference");
    chi2_sig->SetLineColor(kBlue);

    THStack* stack = new THStack();
    stack->Add(chi2_bg);
    stack->Add(chi2_sig);

    new TCanvas();
    stack->Draw("nostack");
    stack->GetXaxis()->SetTitle("#chi^{2}/dof");
    gPad->BuildLegend();
}

void OmegaEtaG::Chi2CutCompare(TTree* tree, double chi2)
{

    new TCanvas();
    TCut cut = Form("EPB_chi2dof<%lf",chi2);

    TH1* all_nocut = Draw(tree, "ggg_IM","",1200,0,1200);
    all_nocut->SetTitle("No Cut");
    all_nocut->SetLineColor(kBlack);

    TH1* bg_nocut = Draw(tree, "ggg_IM", "SigBgFlag==2", 1200,0,1200);
    bg_nocut->SetTitle("Background No Cut");
    bg_nocut->SetLineColor(kGray);

    TH1* all_cut = Draw(tree, "ggg_IM",cut,1200,0,1200);
    all_cut->SetTitle("Cut");
    all_cut->SetLineColor(kRed);

    TH1* bg_cut = Draw(tree, "ggg_IM", cut+"SigBgFlag==2", 1200,0,1200);
    bg_cut->SetTitle("Background Cut");
    bg_cut->SetLineColor(kOrange);

    TH1* sig_cut = Draw(tree, "ggg_IM",cut+"SigBgFlag!=2",1200,0,1200);
    sig_cut->SetTitle("Signal + Reference Cut");
    sig_cut->SetLineColor(kGreen);

    TH1* sig_nocut = Draw(tree, "ggg_IM", "SigBgFlag!=2", 1200,0,1200);
    sig_nocut->SetTitle("Signal + Reference No Cut");
    sig_nocut->SetLineColor(kBlue);

    THStack* stack = new THStack();
    stack->Add(all_nocut);
    stack->Add(all_cut);
    stack->Add(bg_nocut);
    stack->Add(bg_cut);
    stack->Add(sig_nocut);
    stack->Add(sig_cut);
    stack->Draw("nostack");
    gPad->BuildLegend(0.1,0.7,0.5,0.9);
    stack->GetXaxis()->SetTitle("3#gamma IM [MeV]");
}

void OmegaEtaG::gggIM_MM(TTree* tree, const TCut& cut)
{
    new TCanvas();

    TH2* h = Draw(tree, "mmvect_IM:ggg_IM", cut*"TagW", 600, 0, 1200, 750, 500,2000);
    h->SetXTitle("3#gamma IM [MeV]");
    h->SetYTitle("calculated p IM [MeV]");


}

void OmegaEtaG::ggIM(TTree* tree, const TCut& cut)
{
    new TCanvas();

    TH1* h = Draw(tree, "ggIM", cut*"TagW", 1000, 0, 1000);
    h->SetXTitle("2#gamma IM [MeV]");
    h->SetYTitle("#/1MeV");
}

void OmegaEtaG::Analyse(TTree* tree) {

    new TCanvas();

    Chi2Plot(tree);
    Chi2CutCompare(tree,4);

}

void OmegaEtaG::kinfit1() {

    if(gDirectory->cd("OmegaEtaG2") ) {

        TTree* tree = NULL;
        gDirectory->GetObject("tree", tree);

        if(tree)
            Analyse(tree);
    }
}

void OmegaEtaG::Plot(TTree* tree, const double binscale)
{
    const TCut cut_fit("EPB_chi2dof<4 && EPB_status==0");
    const TCut cut_vars("ggg_IM>700&&ggg_IM<840");

    const TCut propmt_randomt("TagW");

    const TCut all = propmt_randomt * (cut_fit + cut_vars);

    const BinSettings IMbins(1600*binscale,   0, 1600);
    const BinSettings MMbins(1600*binscale, 400, 2000);

    const BinSettings MMgggIMbins_X(600*binscale, 0, 1200);
    const BinSettings MMgggIMbins_Y(750*binscale, 500, 2000);


    TH1* ggg_IM_free   = Draw(tree, "ggg_IM",    propmt_randomt,         "3#gamma IM [MeV]", "", IMbins);
    TH1* ggg_IM_fit    = Draw(tree, "ggg_IM",    propmt_randomt*cut_fit, "3#gamma IM [MeV]", "", IMbins);

    TH1* ggIM          = Draw(tree, "ggIM",      all,                    "2#gamma IM [MeV]", "", IMbins);

    TH1* mm_free       = Draw(tree, "mmvect_IM", propmt_randomt,         "MM [MeV]", "", MMbins);
    TH1* mm_fit        = Draw(tree, "mmvect_IM", propmt_randomt*cut_fit, "MM [MeV]", "", MMbins);

    TH2* mm_gggIM_free = Draw(tree, "mmvect_IM:ggg_IM", propmt_randomt,         "3#gamma IM[MeV]", "MM [MeV]", MMgggIMbins_X, MMgggIMbins_Y);
    TH2* mm_gggIM_fit  = Draw(tree, "mmvect_IM:ggg_IM", propmt_randomt*cut_fit, "3#gamma IM[MeV]", "MM [MeV]", MMgggIMbins_X, MMgggIMbins_Y);

    canvas("Cut Varaibles")
            << ggg_IM_free << ggg_IM_fit
            << mm_free << mm_fit
            << drawoption("colz") << mm_gggIM_free << mm_gggIM_fit
            << endc;

    canvas("2g IM")
            << ggIM << endc;
}
