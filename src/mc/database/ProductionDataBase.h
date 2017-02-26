#pragma once

#include <vector>

#include "base/ParticleTypeTree.h"

#include "Math/Interpolator.h"


namespace ant
{
namespace mc
{
namespace data
{

class ProductionDataBase
{
public:
    using XSections_t = std::map<ParticleTypeTreeDatabase::Channel,std::function<double(double)>>;
    struct DataPoint
    {
        const double Energy;
        const double Xsection;
        DataPoint(double energy, double xsection):
            Energy(energy),
            Xsection(xsection){}
    };


    /// Generate a tweaked Interpolator for cross-sections from gichen datapoint
    static std::function<double(double)> MakeInterPolator( const std::vector<DataPoint>& data);

    /// Main Database for cross-sections:
    static const XSections_t XSections;

private:
    /// constructor
    static XSections_t MakeXSections();
};

}
}
}
