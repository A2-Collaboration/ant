#include "OmegaEtaG.h"

#include "TTree.h"
#include "TFile.h"
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
#include "base/std_ext/string.h"

#include "TreeTools.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;
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

    const TCut prompt_random("TagW");

    const TCut all = prompt_random * (cut_fit + cut_vars);

    const BinSettings Ebins (1600*binscale,   0, 1600);

    const BinSettings Chi2bins (250*binscale, 0,   25);
    const BinSettings probbins (250*binscale, 0,   1);

    const BinSettings IMbins(1600*binscale,   0, 1600);
    const BinSettings MMbins(1600*binscale, 400, 2000);

    const BinSettings MMgggIMbins_X(600*binscale, 0, 1200);
    const BinSettings MMgggIMbins_Y(750*binscale, 500, 2000);


    TH1* ggg_IM_free   = Draw(tree, "ggg_IM",    prompt_random,         "3#gamma IM [MeV]", "", IMbins);
    TH1* ggg_IM_fit    = Draw(tree, "ggg_IM",    prompt_random*cut_fit, "3#gamma IM [MeV]", "", IMbins);

    TH1* ggIM          = Draw(tree, "ggIM",      all,                    "2#gamma IM [MeV]", "", IMbins);

    TH1* mm_free       = Draw(tree, "mmvect_IM", prompt_random,         "MM [MeV]", "", MMbins);
    TH1* mm_fit        = Draw(tree, "mmvect_IM", prompt_random*cut_fit, "MM [MeV]", "", MMbins);

    TH2* mm_gggIM_free = Draw(tree, "mmvect_IM:ggg_IM", prompt_random,         "3#gamma IM[MeV]", "MM [MeV]", MMgggIMbins_X, MMgggIMbins_Y);
    TH2* mm_gggIM_fit  = Draw(tree, "mmvect_IM:ggg_IM", prompt_random*cut_fit, "3#gamma IM[MeV]", "MM [MeV]", MMgggIMbins_X, MMgggIMbins_Y);

    auto fit_chi2dof   = Draw(tree, "EPB_chi2dof",     prompt_random,        "#chi^{2}/dof", "", Chi2bins);
    auto fit_prob      = Draw(tree, "EPB_probability", prompt_random,        "probabiliyt", "",  probbins);

//    auto bachelor_free     = Draw(tree, "BachelorE", propmt_randomt,             "E #gamma_{#omega} [MeV]", "", Ebins);
//    auto bachelor_fit      = Draw(tree, "BachelorE", propmt_randomt*cut_fit,     "E #gamma_{#omega} [MeV]", "", Ebins);
    auto bachelor          = Draw(tree, "BachelorE", all,                        "E #gamma_{#omega} [MeV]", "", Ebins);

    map<string, TH1*> pulls;

    for(const auto& particle : {"Proton", "Photon0", "Photon1","Photon2"}) {
        for(const auto& var : {"Phi_pull","Theta_pull","Ek_pull"}) {
            const string branch = formatter() << "EPB_" << particle << "_" << var;
            cout << branch << endl;
            auto h = Draw(tree, branch, prompt_random*cut_fit, "", "", BinSettings(200,-20,20));
            if(h)
                pulls[branch] = h;
        }
    }

    canvas("Kin Fit")
            << fit_chi2dof
            << fit_prob
            << endc;

    canvas("Cut Varaibles")
            << ggg_IM_free << ggg_IM_fit
            << mm_free << mm_fit
            << drawoption("colz") << mm_gggIM_free << mm_gggIM_fit
            << endc;

    canvas("2g IM")
            << ggIM << bachelor << endc;

    canvas cpulls("Pulls");

    for(auto& h : pulls)
        cpulls << h.second;

    cpulls << endc;

}

void OmegaEtaG::DataMC(TFile* mc_file, TFile* data_file, const double mcscale)
{
    canvas c("MC Data");
    for(int i=0;i<8; ++i) {
        TH1* mc   = nullptr;
        TH1* data = nullptr;

        const auto histname = Form("h1d_%d",i);

        mc_file->GetObject(histname, mc);
        if(!mc) {
            cout << "Could not get " << histname << " from mc file" << endl;
            continue;
        }
        cout << "MC: " << mc->GetName() << endl;
        mc->SetLineColor(kBlue);
        const string backup_title(mc->GetTitle());
        mc->SetTitle("mc");
        mc->Scale(mcscale);

        data_file->GetObject(histname, data);
        if(!data) {
            cout << "Could not get " << histname << " from mc file" << endl;
            continue;
        }
        data->SetLineColor(kBlack);
        data->SetTitle("data");


        hstack* stack = new hstack(Form("stack_%d",i), backup_title);
        *stack << mc << data;

        c << drawoption("nostack") << padoption::Legend << *stack;

    }

    c << endc;
}
