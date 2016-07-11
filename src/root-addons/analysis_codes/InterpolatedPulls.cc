#include "InterpolatedPulls.h"

#include "base/std_ext/string.h"
#include <string>
#include "TDirectory.h"
#include "TCanvas.h"
#include <list>
#include <TH2.h>
#include <TH1.h>
#include "TList.h"
#include "TKey.h"
#include <iostream>
#include "base/TH1Ext.h"
#include "THStack.h"
#include "TClass.h"

#include "analysis/plot/root_draw.h"

using namespace std;
using namespace ant;
using namespace ant::std_ext;

template <typename T>
T* getObj(TDirectory* d, const string& name) {
    T* h = nullptr;
    d->GetObject(name.c_str(), h);
    if(!h) {
        cerr << T::Class()->GetName() << " named " << name << " not found!" << endl;
    }
    return h;
}


void InterpolatedPulls::PlotDirectory(TDirectory* dir, const string& prefix, const int cols, const int rows, const string& title, TDirectory* dir2) {

    auto c = new TCanvas();
    c->SetTitle(title.c_str());
    c->Divide(cols, rows, 0, 0);

    for(int y=0; y< rows; ++y) {
        for(int x=0; x< cols; ++x) {

            c->cd( 1 + x + (rows-y-1)*cols);

            const string name = std_ext::formatter() << prefix << "_" << x+1 << "_" << y+1;

            TH1* h1 = nullptr;
            TH1* h2 = nullptr;

            h1 = getObj<TH1>(dir, name);

            if(!h1)
                continue;

            h1->SetStats(false);

            if(dir2) {
                h2 = getObj<TH1>(dir2, name);
                if(!h2)
                    continue;

                h2->SetStats(false);

                THStack* s = new THStack();
                s->Add(h1);
                h1->SetLineColor(kRed);
                s->Add(h2);
                h2->SetLineColor(kBlue);
                s->Draw();

            } else {
                h1->Draw();
            }

        }
    }

}

void InterpolatedPulls::PlotPullSigmas(const string& treename, TDirectory* dir)
{
    canvas c(treename);
    c << drawoption("colz");

    for(const auto& n : {"E","Theta", "Phi"}) {
        TH2* newsigmas = nullptr;
        TH2* oldsigmas = nullptr;
        TH2* pulls = nullptr;

        newsigmas = getObj<TH2>(dir, formatter() << "sigma_" << n);
        oldsigmas = getObj<TH2>(dir, formatter() << "h_sigmas" << n << "_FitZ/h_sigmas" << n << "_z_Mean");
        pulls     = getObj<TH2>(dir, formatter() << "h_pulls" << n <<"_FitZ/h_pulls" << n << "_z_RMS");

        if(!newsigmas || !oldsigmas || !pulls)
            continue;

        TH2Ext::MakeSameZRange({newsigmas, oldsigmas});

        c << newsigmas << oldsigmas << pulls;

    }

    c << endc;

}

void InterpolatedPulls::PlotComparePulls(TDirectory* red, TDirectory* blue)
{
    for(const auto& particle : {"photon", "proton"}) {
        for(const auto& det : {"cb", "taps"}) {
            const string treename = formatter() << "sigma_" << particle << "_" << det;
            TDirectory* d_red = nullptr;
            TDirectory* d_blue = nullptr;

            d_red  = getObj<TDirectory>(red,  treename);
            d_blue = getObj<TDirectory>(blue, treename);

            if(!d_red || !d_blue)
                continue;

            for(const auto& varname : {"E", "Theta", "Phi"}) {
                TDirectory* d_var_red  = nullptr;
                TDirectory* d_var_blue = nullptr;

                const string var_d_name = formatter() << "h_pulls" << varname << "_FitZ";
                const string var_h_name = formatter() << "h_pulls" << varname << "_z";

                d_var_red  = getObj<TDirectory>(d_red,  var_d_name);
                d_var_blue = getObj<TDirectory>(d_blue, var_d_name);

                if(!d_var_red || !d_var_blue)
                    continue;

                PlotDirectory(d_var_red, var_h_name, 15, 10, formatter() << particle << " " << det << " " <<  varname, d_var_blue);
            }
        }
    }
}

void InterpolatedPulls::PlotAllSigams(TDirectory* dir)
{
    for(const auto& particle : {"photon", "proton"}) {
        for(const auto& det : {"cb", "taps"}) {
            const string treename = formatter() << "sigma_" << particle << "_" << det;
            TDirectory* d = nullptr;
            d = getObj<TDirectory>(dir, treename);
            if(d)
                PlotPullSigmas(formatter() << dir->GetName() << ": sigma_" << particle << "_" << det, d);
        }
    }
}
