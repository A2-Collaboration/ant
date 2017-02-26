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

struct ProductionDataBase
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

    static XSections_t MakeXSections();

    static std::function<double(double)> MakeInterPolator( const std::vector<DataPoint>& data);

    static const XSections_t XSections;

};

}
}
}
