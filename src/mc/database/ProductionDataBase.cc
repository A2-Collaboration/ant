#include "ProductionDataBase.h"

#include "base/interval.h"

#include "detail/gp_pPi0.h"
#include "detail/compton.h"


using namespace std;
using namespace ant::mc::data;


const ProductionDataBase::XSections_t ProductionDataBase::XSections = MakeXSections();


ProductionDataBase::XSections_t ProductionDataBase::MakeXSections()
{
    XSections_t db;

    db.insert(gp_ppi0);
    db.insert(compton);



    // Estimations for EPT-Range only !!!!
    auto makeBox = [](const IntervalD& range, double xsecEstimate){
        return [range,xsecEstimate](double Egamma){
            return range.Contains(Egamma) ?  xsecEstimate : 0.0;
        };
    };
    db.insert({ParticleTypeTreeDatabase::Channel::ThreePi0_6g,
               makeBox({1410,1600},1.2) });
    db.insert({ParticleTypeTreeDatabase::Channel::gp_pPiPPiMPi0,
               makeBox({1410,1600},25.0) });
    db.insert({ParticleTypeTreeDatabase::Channel::gp_pPiPPiMPi0Pi0,
               makeBox({1410,1600},18.0) });

    return db;
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

