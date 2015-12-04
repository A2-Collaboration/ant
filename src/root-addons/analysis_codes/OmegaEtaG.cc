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

#include "TreeTools.h"

using namespace std;

void OmegaEtaG::Chi2Plot(TTree* tree)
{
    new TCanvas("chi2_plot","Kinfit Chi2");
    TH1* chi2_bg = Draw(tree, "kinfit_chi2","rf==2", 1000,0,10);
    chi2_bg->SetTitle("#chi^{2} Background");
    chi2_bg->SetLineColor(kRed);

    TH1* chi2_sig = Draw(tree, "kinfit_chi2","rf!=2", 1000,0,10);
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
    TCut cut = Form("kinfit_chi2<%lf",chi2);

    TH1* all_nocut = Draw(tree, "ggg_IM","",1200,0,1200);
    all_nocut->SetTitle("No Cut");
    all_nocut->SetLineColor(kBlack);

    TH1* bg_nocut = Draw(tree, "ggg_IM", "rf==2", 1200,0,1200);
    bg_nocut->SetTitle("Background No Cut");
    bg_nocut->SetLineColor(kGray);

    TH1* all_cut = Draw(tree, "ggg_IM",cut,1200,0,1200);
    all_cut->SetTitle("Cut");
    all_cut->SetLineColor(kRed);

    TH1* bg_cut = Draw(tree, "ggg_IM", cut+"rf==2", 1200,0,1200);
    bg_cut->SetTitle("Background Cut");
    bg_cut->SetLineColor(kOrange);

    TH1* sig_cut = Draw(tree, "ggg_IM",cut+"rf!=2",1200,0,1200);
    sig_cut->SetTitle("Signal + Reference Cut");
    sig_cut->SetLineColor(kGreen);

    TH1* sig_nocut = Draw(tree, "ggg_IM", "rf!=2", 1200,0,1200);
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

    TH2* h = Draw(tree, "calcp_IM:ggg_IM", cut, 600, 0, 1200, 750, 500,2000);
    h->SetXTitle("3#gamma IM [MeV]");
    h->SetYTitle("calculated p IM [MeV]");


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
