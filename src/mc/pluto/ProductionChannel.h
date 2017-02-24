#pragma once

#include <vector>

#include "base/ParticleTypeTree.h"

#include "Math/Interpolator.h"


namespace ant
{
namespace mc
{
namespace pluto
{

struct ChannelDataBase
{
    using XSections_t = std::map<ParticleTypeTreeDatabase::Channel,std::function<double(double)>>;

    struct DataPoint
    {
        const double Energy;
        const double Xsection;
        DataPoint(double energy, double xsection):
            Energy(energy),
            Xsection(xsection){}
    };

    static std::function<double(double)> MakeInterPolator( const std::vector<DataPoint>& data)
    {
        std::vector<double> dataE;
        std::vector<double> dataXsec;

        for (const auto& d: data)
        {
            dataE.emplace_back(d.Energy);
            dataXsec.emplace_back(d.Xsection);
        }
        return [dataE,dataXsec] (double energy)
        {
            return ROOT::Math::Interpolator(dataE, dataXsec).Eval(energy);
        };
    }

    static XSections_t MakeXSections();
    static const XSections_t XSections;

};

}
}
}
