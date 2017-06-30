#include "pi0.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include "root-addons/analysis_codes/Math.h"
#include "TAxis.h"
#include "analysis/plot/RootDraw.h"
#include "TGraph.h"
#include "base/ParticleType.h"
#include "analysis/plot/HistogramFactory.h"
#include "base/math_functions/CrystalBall.h"
#include "base/math_functions/AsymGaus.h"
#include "base/std_ext/math.h"
#include "base/std_ext/string.h"

#include "base/std_ext/memory.h"
#include "analysis/utils/MCWeighting.h"
#include "root-addons/analysis_codes/hstack.h"
#include "TNtupleD.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "base/TH_ext.h"
#include "base/PhysicsMath.h"

#include "TROOT.h"
#include "TFile.h"

using namespace std;
using namespace ant;



void ant::Pi0::plotSigmaTheta(const bool showRel)
{
    auto canvas = new TCanvas();
    TH1D* hstd = nullptr;



    auto files = gROOT->GetListOfFiles();
    if (files->GetEntries() == 0)
    {
        cout << "No files present!" << endl;
        return;
    }

    vector<TH1D*> hists(files->GetEntries());

    for (int i = 0 ; i < files->GetEntries() ; ++i)
    {
        auto file = dynamic_cast<TFile*>(files->At(i));
        if (!file) continue;

        TH1D* hist  = nullptr;

        file->GetObject("SigmaTheta",hist);
        if (!hist) continue;
        hists.at(i) = hist;

        const string name  = std_ext::basename(file->GetName());

        hist->SetName(name.c_str());

        if (std_ext::contains(name,"dek<20_emb0.05_Pi0PIDVeto==0.root"))
        {
            hstd = hist;
            continue;
        }

        if (std_ext::contains(name,"dek==0"))
        {
            hist->SetLineColor(kRed);
        }
        if (std_ext::contains(name,"dek<20"))
        {
            hist->SetLineColor(kBlue);
        }

        if (std_ext::contains(name,"dek<50"))
        {
            hist->SetLineColor(kBlack);
        }

        auto ya = hist->GetYaxis();
        ya->SetRangeUser(0,2.3);
        ya->SetNdivisions(5,5,0,true);
        ya->SetTitleOffset(0.77);
        ya->SetTitleSize(0.05);
        ya->SetLabelSize(0.05);
        auto xa = hist->GetXaxis();
        xa->SetNdivisions(5,5,0,true);
        xa->SetTitleOffset(0.65);
        xa->SetTitleSize(0.05);
        xa->SetLabelSize(0.05);

    }


    canvas->SetCanvasSize(600 + showRel * 700 ,600);
    if (showRel) canvas->Divide(2);


    auto pad1 = canvas->cd(1);
    pad1->SetGridx();
    pad1->SetGridy();

    for (auto h: hists)
    {
        if (!(h == hstd))
            h->Draw("same");
    }
    if (hstd)
    {
        hstd->SetLineWidth(3);
        hstd->SetLineColor(kGreen);
        hstd->Draw("same");
    }




    if (showRel)
    {
        TH1D* clone = dynamic_cast<TH1D*>(hstd->Clone("clone"));
        vector<TH1D*> hists_rel(hists.size());
        transform(hists.begin(),hists.end(),hists_rel.begin(),
                  [](TH1D* h)
                  {
                      const string name = std_ext::formatter() << "rel_" << h->GetName();
                      auto clone = dynamic_cast<TH1D*>(h->Clone(name.c_str()));
                      clone->SetYTitle("relative change");
                      auto ya = clone->GetYaxis();
                      ya->SetRangeUser(0.7,1.4);
                      ya->SetNdivisions(5,5,0,true);
                      ya->SetTitleOffset(1.00);
                      ya->SetTitleSize(0.05);
                      ya->SetLabelSize(0.05);

                      return clone;
                  } );

        auto pad2 = canvas->cd(2);
        pad2->SetGridx();
        pad2->SetGridy();
        bool same  = false;
        for (auto h: hists_rel)
        {
            const string drwo = same ? "same" : "";
            same = true;

            if (clone)
                h->Divide(clone);
            h->Draw("same");
        }
    }

}

void Pi0::plotComparison()
{
}
