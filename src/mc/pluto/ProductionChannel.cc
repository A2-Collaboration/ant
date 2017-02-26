#include "ProductionChannel.h"

#include "database/gp_ppi0.h"


using namespace std;
using namespace ant::mc::pluto;


const ProductionDataBase::XSections_t ProductionDataBase::XSections = MakeXSections();


ProductionDataBase::XSections_t ProductionDataBase::MakeXSections()
{
    XSections_t x;

    x.insert(gp_ppi0);

    return x;
}



std::function<double (double)> ProductionDataBase::MakeInterPolator(const std::vector<ProductionDataBase::DataPoint>& data)
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
        // don't crash for unknown energies, just use last known value
        if (energy < dataE.front()){
            return  dataXsec.front();
        }
        if (energy > dataE.back()){
            return  dataXsec.back();
        }

        auto xsec = ROOT::Math::Interpolator(dataE, dataXsec).Eval(energy);
        //quickfix for Interpolator smoothing into negative numbers
        if (xsec < 0 ) xsec = 0;

        return xsec;
    };
}

