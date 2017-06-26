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
#include "base/PlotExt.h"
#include "THStack.h"
#include "TClass.h"
#include <list>
#include "tree/TCandidate.h"
#include "tree/TParticle.h"
#include "Rtypes.h"
#include "TROOT.h"
#include "analysis/utils/uncertainties/Interpolated.h"
#include "analysis/utils/uncertainties/Optimized.h"
#include "TGraph.h"
#include "analysis/plot/RootDraw.h"
#include "TLine.h"

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

void InterpolatedPulls::PlotDirectory(list<TDirectory*> dirs, const string& prefix, const string& title) {

    const vector<Color_t> colors = {kRed, kBlue, kGreen, kBlack};

    if(dirs.empty())
        return;

    auto hist2d = getObj<TH2>(dirs.front(), std_ext::formatter() << prefix << "_Mean");
    const auto cols = hist2d->GetNbinsX();
    const auto rows = hist2d->GetNbinsY();

    auto c = new TCanvas();
    c->SetTitle(title.c_str());
    c->Divide(cols, rows, 0, 0);

    for(int y=0; y< rows; ++y) {
        for(int x=0; x< cols; ++x) {

            c->cd( 1 + x + (rows-y-1)*cols);

            const string name = std_ext::formatter() << prefix << "_" << x+1 << "_" << y+1;

            THStack* s = new THStack();

            size_t c=0;
            for(auto d : dirs) {

                auto h1 = getObj<TH1>(d, name);
                if(h1) {

                    s->Add(h1);
                    h1->SetLineColor(colors.at(c));
                    h1->SetStats(false);
                    c = (c+1) % colors.size();
                }
            }

            s->Draw("nostack");


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

                PlotDirectory( {d_var_red, d_var_blue}, var_h_name, formatter() << particle << " " << det << " " <<  varname);
            }
        }
    }
}

void InterpolatedPulls::PlotComparePulls2()
{
    auto files = gROOT->GetListOfFiles();

    for(const auto& particle : {"photon", "proton"}) {
        for(const auto& det : {"cb", "taps"}) {
            const string treename = formatter() << "sigma_" << particle << "_" << det;

            list<TDirectory*> dirs;
            for(int i=0; i<files->GetEntries(); ++i) {
                auto f = dynamic_cast<TDirectory*>(files->At(i));
                if(!f)
                    throw std::runtime_error("Found something odd in list of files");

                auto d = getObj<TDirectory>(f,  treename);

                if(d)
                    dirs.push_back(d);
            }

            if(dirs.empty())
                continue;

            for(const auto& varname : {"E", "Theta", "Phi"}) {

                list<TDirectory*> var_dirs;


                const string var_d_name = formatter() << "h_pulls" << varname << "_FitZ";
                const string var_h_name = formatter() << "h_pulls" << varname << "_z";

                for( auto d: dirs) {
                    auto vd = getObj<TDirectory>(d,  var_d_name);
                    if(vd)
                        var_dirs.push_back(vd);

                }

                PlotDirectory(var_dirs, var_h_name, formatter() << particle << " " << det << " " <<  varname);
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

    const double xmin = det == Detector_t::Type_t::TAPS ? 0.9 : -0.9;
    const double xmax = det == Detector_t::Type_t::TAPS ? 1.0 :  0.9;

    const string E_title = formatter() << particle.Name() << " " << Detector_t::ToString(det) << " Ek, interpolated";
    TH2D* h_Ek     = new TH2D("", E_title.c_str(),     100, xmin, xmax, 100, 0, 1000);

    const string Theta_title = formatter() << particle.Name() << " " << Detector_t::ToString(det) << " Theta, interpolated";
    TH2D* h_Theta = new TH2D("",Theta_title.c_str(),  100, xmin, xmax, 100, 0, 1000);

    const string Phi_title = formatter() << particle.Name() << " " << Detector_t::ToString(det) << " Phi, interpolated";
    TH2D* h_Phi   = new TH2D("", Phi_title.c_str(),   100, xmin, xmax, 100, 0, 1000);

    for(int x=1; x<h_Ek->GetNbinsX(); ++x) {
        for(int y=1; y<h_Ek->GetNbinsY(); ++y) {
            const auto Theta = acos(h_Ek->GetXaxis()->GetBinCenter(x));
            const auto E     = h_Ek->GetYaxis()->GetBinCenter(y);

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

            h_Ek->SetBinContent(     x, y, v.sigmaEk);
            h_Theta->SetBinContent( x, y, v.sigmaTheta);
            h_Phi->SetBinContent(   x, y, v.sigmaPhi);
        }
    }

    canvas() << drawoption("colz") << h_Ek << h_Theta << h_Phi << endr << endc;
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

inline TGraph* makeGraph(const std::vector<vec2>& points) {
    if(points.empty())
        return nullptr;
    auto g = new TGraph(int(points.size()));
    for(size_t p=0; p<points.size(); ++p) {
        const auto& v = points.at(p);
        g->SetPoint(int(p), v.x, v.y);
    }
    return g;
}

TCanvas* InterpolatedPulls::ConvergencePlots(std::list<TH2*> hists)
{
    if(hists.empty())
        return nullptr;

    const auto cols = hists.front()->GetNbinsX();
    const auto rows = hists.front()->GetNbinsY();

    TCanvas* c = new TCanvas();
    c->Divide(cols, rows, 0, 0);

    for(int x=0; x<cols; ++x) {
        for(int y=0; y<rows; ++y) {
            c->cd( 1 + x + (rows-y-1)*cols);

            vector<vec2> pnts;

            int p=0;
            for(auto h : hists) {
                const auto v = h->GetBinContent(x+1,y+1);
                if(isfinite(v))
                    pnts.emplace_back(vec2(p,v));
                ++p;
            }
            auto g = makeGraph(pnts);
            if(g) {
                g->GetYaxis()->SetRangeUser(0.5, 1.5);
                g->SetTitle("");
                g->Draw("ALP");
                TLine* l = new TLine(g->GetXaxis()->GetXmin(), 1, g->GetXaxis()->GetXmax(), 1);
                l->SetLineColor(kGray);
                l->Draw();
            }
        }
    }

    return c;

}

void InterpolatedPulls::ConvergencePlots2()
{
    auto files = gROOT->GetListOfFiles();

    for(const auto& particle : {"photon", "proton"}) {
        for(const auto& det : {"cb", "taps"}) {
            const string treename = formatter() << "sigma_" << particle << "_" << det;

            list<TDirectory*> dirs;
            for(int i=0; i<files->GetEntries(); ++i) {
                auto f = dynamic_cast<TDirectory*>(files->At(i));
                if(!f)
                    throw std::runtime_error("Found something odd in list of files");

                auto d = getObj<TDirectory>(f,  treename);

                if(d)
                    dirs.push_back(d);
            }

            if(dirs.empty())
                continue;

            for(const auto& varname : {"E", "Theta", "Phi"}) {

                list<TH2*> var_hists;


                const string var_d_name = formatter() << "h_pulls" << varname <<"_FitZ/h_pulls" << varname << "_z_RMS";

                for( auto d: dirs) {
                    auto vd = getObj<TH2>(d,  var_d_name);
                    if(vd)
                        var_hists.push_back(vd);

                }

                auto c = ConvergencePlots(var_hists);
                const string title = formatter() << particle << " " << det << " " << varname;
                c->SetTitle(title.c_str());
            }
        }
    }
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

void ConvergencePlot::Plot(const double min, const double max, const int ww, const int wh)
{

    for(const auto& particle : {"photon","proton"}) {
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

                if(min != max) {
                    TH2Ext::MakeSameZRange(hists, {min,max});
                }

                TH2Ext::MakeSameZRange(hists);

                const string title = formatter() << particle << " " << det << " " << n;

                TCanvas* c = new TCanvas();
                c->SetTitle(title.c_str());
                c->Divide(int(hists.size()),1);

                if(ww>0 && wh>0) {
                    c->SetCanvasSize(ww,wh);
                }

                int p=1;
                for(auto h : hists) {
                    c->cd(p++);
                    h->Draw("colz");
                }

                const string fname = formatter() << "conv_col_" << particle << "_" << det << "_" << n;
                auto save_multiimages = [] (TCanvas* c, const char* basename) {
                    c->SaveAs(Form("%s.pdf", basename));
                    c->SaveAs(Form("%s.png", basename));
                    c->SaveAs(Form("%s.root", basename));
                };
                save_multiimages(c, fname.c_str());

            }
        }
    }
}


