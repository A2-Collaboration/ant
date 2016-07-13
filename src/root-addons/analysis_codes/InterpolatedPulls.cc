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
#include <list>
#include "tree/TCandidate.h"
#include "tree/TParticle.h"
#include "analysis/plot/HistogramFactories.h"

#include "analysis/utils/Uncertainties.h"

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

std::pair<int,int> findMaxXY(TDirectory* dir, const string& prefix) {
    int nx=-1;
    int x=-1;
    bool row_ok=false;
    bool cell_ok=false;
    int y=-1;

    do {
        ++y;
        row_ok=false;
        x=0;
        do {

            cell_ok = getObj<TH1>(dir, std_ext::formatter() << prefix << "_" << x +1 << "_" << y+1) != nullptr;

            if(cell_ok) {
                row_ok=true;
                ++x;
            }

        } while (cell_ok);

        if(nx==-1) {
            nx=x;
        } else if(row_ok) {
            if(nx != x) {
                cerr << "Row length mismatch! previous had " << nx << ", this has " << x << endl;
                return {-1,-1};
            }
        }


    } while (row_ok);

    return {nx, y};

}

void InterpolatedPulls::PlotDirectory(TDirectory* dir, const string& prefix, const string& title, TDirectory* dir2) {

    auto hist2d = getObj<TH2>(dir, std_ext::formatter() << prefix << "_Mean");
    const auto cols = hist2d->GetNbinsX();
    const auto rows = hist2d->GetNbinsY();

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

                PlotDirectory(d_var_red, var_h_name, formatter() << particle << " " << det << " " <<  varname, d_var_blue);
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

void scanModel(analysis::utils::UncertaintyModelPtr model,  const Detector_t::Type_t det, const ParticleTypeDatabase::Type& particle) {

    const double xmin = det & Detector_t::Type_t::TAPS ? 0.9 : -0.9;
    const double xmax = det & Detector_t::Type_t::TAPS ? 1.0 :  0.9;

    const string E_title = formatter() << particle.Name() << " " << Detector_t::ToString(det) << " E, interpolated";
    TH2D* h_E     = new TH2D("", E_title.c_str(),     100, xmin, xmax, 100, 0, 1000);

    const string Theta_title = formatter() << particle.Name() << " " << Detector_t::ToString(det) << " Theta, interpolated";
    TH2D* h_Theta = new TH2D("",Theta_title.c_str(),  100, xmin, xmax, 100, 0, 1000);

    const string Phi_title = formatter() << particle.Name() << " " << Detector_t::ToString(det) << " Phi, interpolated";
    TH2D* h_Phi   = new TH2D("", Phi_title.c_str(),   100, xmin, xmax, 100, 0, 1000);

    for(int x=1; x<h_E->GetNbinsX(); ++x) {
        for(int y=1; y<h_E->GetNbinsY(); ++y) {
            const auto Theta = acos(h_E->GetXaxis()->GetBinCenter(x));
            const auto E     = h_E->GetYaxis()->GetBinCenter(y);

            auto cand = make_shared<TCandidate>(det,
                                                E,
                                                Theta,
                                                0.0, // Phi
                                                0.0, // Time
                                                1,   // cluster size
                                                0.0, // veto
                                                0.0, // tracker
                                                TClusterList{}   // Cluster List
                                                );

            const auto part = TParticle(particle, cand);

            const auto v = model->GetSigmas(part);

            h_E->SetBinContent(     x, y, v.sigmaE);
            h_Theta->SetBinContent( x, y, v.sigmaTheta);
            h_Phi->SetBinContent(   x, y, v.sigmaPhi);
        }
    }

    canvas() << drawoption("colz") << h_E << h_Theta << h_Phi << endr << endc;
}

void InterpolatedPulls::TestInterpolation(const string& filename)
{
    auto dflt = make_shared<analysis::utils::UncertaintyModels::Optimized_Oli1>();
    auto model = make_shared<analysis::utils::UncertaintyModels::Interpolated>(dflt);
    model->LoadSigmas(filename);

    scanModel(model, Detector_t::Type_t::CB,   ParticleTypeDatabase::Photon);
    scanModel(model, Detector_t::Type_t::TAPS, ParticleTypeDatabase::Photon);
    scanModel(model, Detector_t::Type_t::CB,   ParticleTypeDatabase::Proton);
    scanModel(model, Detector_t::Type_t::TAPS, ParticleTypeDatabase::Proton);

}

class DirList : public std::list<TDirectory*> {
    using std::list<TDirectory*>::list;
};

ConvergencePlot::ConvergencePlot():
    list(new DirList())
{}

ConvergencePlot::~ConvergencePlot()
{
    delete list;
}

void ConvergencePlot::Add(TDirectory* dir)
{
    list->push_back(dir);
}

void ConvergencePlot::Plot()
{

    for(const auto& particle : {"photon"}) {
        for(const auto& det : {"cb", "taps"}) {
            for(const auto& n : {"E","Theta", "Phi"}) {

                vector<TH2*> hists;

                for(auto dir : *list) {

                    TH2* pulls = nullptr;

                    pulls     = getObj<TH2>(dir, formatter() << "sigma_" << particle << "_" << det << "/h_pulls" << n <<"_FitZ/h_pulls" << n << "_z_RMS");

                    if(!pulls)
                        return;

                    hists.push_back(pulls);


                }

                TH2Ext::MakeSameZRange(hists);

                const string title = formatter() << particle << " " << det << " " << n;

                TCanvas* c = new TCanvas();
                c->SetTitle(title.c_str());
                c->Divide(hists.size(),1);

                int p=1;
                for(auto h : hists) {
                    c->cd(p++);
                    h->Draw("colz");
                }

            }
        }
    }
}


