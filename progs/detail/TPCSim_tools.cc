#include "TPCSim_tools.h"
#include "TRandom2.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2D.h"
#include "base/PlotExt.h"
#include "base/std_ext/math.h"

#include <algorithm>
#include <iostream>

using namespace ant;
using namespace std;

namespace TPCSim {

vector<ant::vec2> generatePoints(const double z0, const double theta,
                                         const resolution_t &prop,
                                         const tpcproperties& tpc)
{
    const track_t t({0,z0},{cos(std_ext::degree_to_radian(90.0)-theta),sin(std_ext::degree_to_radian(90.0)-theta)});

    vector<vec2> points;
    points.reserve(size_t(tpc.nRings));


    for(int ring = 0; ring < tpc.nRings; ++ring) {
        const double r = tpc.rin + (ring+0.5)*tpc.ringWidth();
        const auto point = generatePoint(r,t,prop);

        if(tpc.length.Contains(point.y))
            points.emplace_back(point);
    }

    return points;
}

vec2 generatePoint(const double r, const track_t &track, const resolution_t &prop)
{
    const auto z = track(r);
    return {r, gRandom->Gaus(z, prop.dl)};
}

TH2D*draw::makeCanvas(const tpcproperties& tpc)
{
    const auto range = tpc.CBradius + tpc.CBradius / 10.;
    auto canvas = new TH2D("","", 100, -range, range,
                                  100, -range, range);
    canvas->SetStats(false);
    canvas->GetXaxis()->SetTitle("r [cm]");
//    canvas->GetXaxis()->SetRangeUser(-3,15); // doesn't do anything
    canvas->GetYaxis()->SetTitle("z [cm]");
    return canvas;
}

TGraphErrors* draw::makeGraph(const std::vector<vec2>& points, const resolution_t &res, const tpcproperties &tpc)
{
    auto g = new TGraphErrors(int(points.size()));

    for(const auto& point : points) {
        const auto error = getUncertainties(point, res, tpc);
        GraphExt::FillGraphErrors(g, point.x, point.y, error.x, error.y);
    }
    return g;
}
TF1* draw::makeFitTF1(const trackFitter_t& tFitter)
{
//    const auto npts = int(tFitter.Fitted_Rs.size());
//    auto g = new TGraphErrors(npts);

//    for(auto i = 0 ; i < npts ; i++) {
//        GraphExt::FillGraph(g,tFitter.Fitted_Rs.at(i).Value(),
//                              tFitter.Fitted_Zs.at(i).Value());
//    }

    auto fitfkt = [&tFitter](double* r, double*)
    {
        return tFitter.a.Value() + tFitter.b.Value() * r[0];
    };

    auto f = new TF1("track",fitfkt,0,100,0);

    return f;
}

std::list<TGraph*> draw::makeScene(const tpcproperties& tpc)
{
    list<TGraph*> list;
    list.push_back(tpc.getOutline());

    //Target
    list.push_back(
                [] () {
        auto g = new TGraph(5);
        constexpr auto l=10.0;
        constexpr auto d=2.0;
        g->SetPoint(0,d/2,l/2);
        g->SetPoint(1,-d/2,l/2);
        g->SetPoint(2,-d/2,-l/2);
        g->SetPoint(3,d/2,-l/2);
        g->SetPoint(4,d/2,l/2);
        return g;
    }());

    //CB
    list.push_back(
                [] () {
                        constexpr auto np = 180;
                        constexpr auto r = 24.5;
                        auto g = new TGraph(np+1);
                        for(int i=0; i<np; ++i) {
                            const auto phi = std_ext::degree_to_radian(360.0/np*i);
                            g->SetPoint(i,r*cos(phi),r*sin(phi));
                        }
                        g->SetPoint(np, r, 0.0);
                        return g;
                    }());
    return list;
}


vec2 getUncertainties(const vec2&, const resolution_t &res, const tpcproperties &tpc)
{
    //        row (fixed)    , longitudinal -> sigma_z
    return {tpc.ringWidth()/2, res.dl};
}

ant::interval<double> tpcproperties::tpcInCB(const double l, const double rout) {
    const auto d = sqrt(CBradius*CBradius - rout*rout);
    return {d-l,d};
}

trackFitter_t::trackFitter_t(const vector<Value_t>& points_r,
                             const vector<Value_t>& points_z):
    Fitted_Rs(points_r), Fitted_Zs(points_z)
{
    APLCON::Fitter<Value_t, Value_t,std::vector<Value_t>,std::vector<Value_t>> fitter;

    auto residuals = [] (const Value_t& a, const Value_t& b, const vector<Value_t>& r, const vector<Value_t>& z) {
        vector<double> residuals(z.size());
        transform(r.begin(), r.end(), z.begin(), residuals.begin(),
                  [&a, &b] (const double& r_i, const double& z_i) {
            return a + b*r_i - z_i;
        });
        return residuals;
    };

    Result = fitter.DoFit(a, b, Fitted_Rs, Fitted_Zs, residuals);
}


TGraph *tpcproperties::getOutline() const
{
    auto g = new TGraph(5);
    g->SetPoint(0,rin,length.Start());
    g->SetPoint(1,rout,length.Start());
    g->SetPoint(2,rout,length.Stop());
    g->SetPoint(3,rin,length.Stop());
    g->SetPoint(4,rin,length.Start());
    return g;
}



tpcproperties::tpcproperties(const double len):
    length(tpcInCB(len,rout))
{}

}

