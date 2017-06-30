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

void ant::Pi0::plotSigmaTheta(const bool rel)
{
    auto canvas = new TCanvas();
    TH1D* hstd = nullptr;

    canvas->SetCanvasSize(600,600);
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
        }

    }



    TH1D* clone = dynamic_cast<TH1D*>(hstd->Clone("clone"));

    for (auto h: hists)
    {
        if (rel && clone)
            h->Divide(clone);
        if (h == hstd)
            continue;
        h->Draw("same");
    }
    if (hstd)
    {
        hstd->SetLineColor(kRed);
        hstd->SetLineWidth(3);
        hstd->Draw("same");
    }

}
