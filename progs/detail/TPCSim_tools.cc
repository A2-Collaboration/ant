#include "TPCSim_tools.h"
#include "TRandom2.h"
#include "TGraphErrors.h"
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

TGraphErrors* makeGraph(const std::vector<vec2>& points, const resolution_t &res, const tpcproperties &tpc)
{
    auto g = new TGraphErrors(int(points.size()));

    for(const auto& point : points) {
        const auto error = getUncertainties(point, res, tpc);
        GraphExt::FillGraphErrors(g, point.x, point.y, error.x, error.y);
    }
    return g;
}
TGraph* makeGraphFitted(const trackFitter_t& tFitter)
{
    const auto npts = int(tFitter.Fitted_Rs.size());
    auto g = new TGraphErrors(npts);

    for(auto i = 0 ; i < npts ; i++) {
        GraphExt::FillGraph(g,tFitter.Fitted_Rs.at(i).Value(),
                            tFitter.Fitted_Zs.at(i).Value());
    }
    return g;
}


vec2 getUncertainties(const vec2&, const resolution_t &res, const tpcproperties &tpc)
{
    return {tpc.ringWidth()/2, res.dt};
}

ant::interval<double> tpcproperties::tpcInCB(const double l, const double rout) {
    const auto d = sqrt(CBradius*CBradius - rout-rout);
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

