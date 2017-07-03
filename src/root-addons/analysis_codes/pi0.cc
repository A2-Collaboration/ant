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

using data_t = vector<pair<double,double>>;
auto makeGraph = [] (const data_t& data,const string& name, const Color_t color, const Style_t style)
{
    auto g = new TGraph();
    g->SetName(name.c_str());
    g->SetMarkerColor(color);
    g->SetMarkerStyle(style);
    g->SetFillColor(kWhite);
    g->SetLineColor(kWhite);

    for (const auto& p: data)
        g->SetPoint(g->GetN(),p.first,p.second);
    return g;
};

const data_t cbElsaTAPS2011_1450_1500 =
{
    { -0.8677419354838705, 1.7999999999999998 },
    { -0.7718963831867054, 1.5696969696969694 },
    { -0.6687683284457471, 1.1454545454545455 },
    { -0.5810850439882698, 0.8909090909090907 },
    { -0.4930107526881722, 0.733333333333333  },
    { -0.404496578690126, 0.684848484848485   },
    { -0.29946236559139816, 0.733333333333333 },
    { 0.4487292277614854, 0.2848484848484847  },
    { 0.5467253176930598, 0.5878787878787874  },
    { 0.6530303030303024, 0.9515151515151516  },
    { 0.7601173020527847, 1.5090909090909088  },
    { 0.8517106549364613, 2.224242424242424   },
    { 0.9571358748778103, 2.369696969696969   }
};
const data_t cbElsa2005_1450_1500 =
{
    { -0.75, 0.8363095238095353                  },
    { -0.6547619047619069, 0.6024659863945696    },
    { -0.5476190476190492, 0.3326955782313048    },
    { -0.4523809523809543, 0.42028061224490987   },
    { -0.34523809523809845, 0.4897959183673586   },
    { -0.2500000000000018, 0.6130952380952501    },
    { -0.15476190476190688, 0.6471088435374273   },
    { -0.047619047619051, 0.8237670068027332     },
    { 0.0476190476190439, 0.5720663265306243     },
    { 0.15476190476190155, 0.4987244897959302    },
    { 0.24999999999999645, 0.264880952380965     },
    { 0.3452380952380949, 0.04889455782314167    },
    { 0.4404761904761898, -0.006377551020395611  },
    { 0.5357142857142847, 0.06335034013606711    },
    { 0.6428571428571423, 0.38286564625851627    },
    { 0.7380952380952372, 0.791879251700693      },
    { 0.8452380952380949, 1.8435374149659995     }

};
const data_t sergey =
{
   { -0.8342828111205502, 1.1220328819149465   },
   { -0.7657621641430545, 1.0011868472340213   },
   { -0.7016599203771258, 0.7898593273699466   },
   { -0.6353588859460628, 0.6306281545150165   },
   { -0.5623615344893641, 0.5618722202788464   },
   { -0.5008203821231654, 0.5370233871370744   },
   { -0.43707211033931326, 0.5587864174703912  },
   { -0.36652757649959167, 0.6052108840137922  },
   { -0.3028709209934535, 0.6873032332217279   },
   { -0.23468342411673548, 0.7858365399030531  },
   { -0.17100178235303898, 0.8514754385088197  },
   { -0.1027726417137238, 0.922586327519864    },
   { -0.04131894124897961, 0.9553245714856826  },
   { 0.03161594456382355, 0.9277022637549346   },
   { 0.10008661902620242, 0.8397631302783468   },
   { 0.1640389452467812, 0.727156314027285     },
   { 0.23027334965768853, 0.611801009444805    },
   { 0.2965952059700503, 0.4388586277547346    },
   { 0.3628087884996587, 0.3372145320073956    },
   { 0.4289224259990341, 0.30138423866873154   },
   { 0.4903677977112588, 0.33960696616860697   },
   { 0.5631361084736732, 0.4216743291189844    },
   { 0.6243191244815356, 0.6326582879416318    },
   { 0.6897664617793551, 1.0355866773274698    },
   { 0.7596905035563775, 1.490605167158063     },
   { 0.827253343994137, 2.0004747388936086     }
};

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
    auto canvas = new TCanvas();

    canvas->SetCanvasSize(600,600);


    auto files = gROOT->GetListOfFiles();
    if (files->GetEntries() == 0)
    {
        cout << "No files present!" << endl;
        return;
    }
    vector<TH1D*> hists(files->GetEntries());


    bool first = true;
    for (int i = 0 ; i < files->GetEntries() ; ++i)
    {
        auto file = dynamic_cast<TFile*>(files->At(i));
        if (!file) continue;

        TH1D* hist  = nullptr;

        file->GetObject("SigmaTheta",hist);
        if (!hist) continue;
        hists.at(i) = hist;

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

        hist->Draw("same");
        if ( first )
        {
            makeGraph(cbElsa2005_1450_1500,"CB-Elsa 2005",kRed,kCircle)->Draw("p same");
            makeGraph(cbElsaTAPS2011_1450_1500,"CB-Elsa TAPS 2011",kBlack,25)->Draw("p same");
            makeGraph(sergey,"CB MAMI 2015",kGreen,kStar)->Draw("p same");
            canvas->BuildLegend();
            first = false;
        }
    }






}
