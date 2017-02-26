#include "ProductionDataBase.h"

#include "base/interval.h"
#include "base/Logger.h"


#include "detail/compton.h"

#include "detail/gp_pPi0.h"
#include "detail/gp_pPi0Pi0.h"

#include "detail/gp_pEta.h"
#include "detail/gp_pEtaPi0.h"

#include "detail/gp_pEtaPrime.h"

#include "detail/gp_pOmega.h"
#include "detail/gp_pRho.h"

#include "detail/gp_nPiP.h"


using namespace std;
using namespace ant::mc::data;


const ProductionDataBase::XSections_t ProductionDataBase::XSections = MakeXSections();


ProductionDataBase::XSections_t ProductionDataBase::MakeXSections()
{
    // helper for constant cross-sections:
    auto makeBox = [](const IntervalD& range, double xsecEstimate){
        return [range,xsecEstimate](double Egamma){
            return range.Contains(Egamma) ?  xsecEstimate : 0.0;
        };
    };

    XSections_t db;

    // insert new channels here:
    for (const auto& ch: {
         compton,
         gp_pPi0, gp_pPi0Pi0,
         gp_pEta, gp_pEtaPi0,
         gp_pEtaPrime,
         gp_pOmega, gp_pRho,
         gp_nPiP}
         )
    {
        db.insert(ch);
    }

    // Estimations for EPT-Range only !!!!
    IntervalD EPTrange({1410,1600});
    db.insert({ParticleTypeTreeDatabase::Channel::gp_p3Pi0,
               makeBox(EPTrange,1.2) });
    db.insert({ParticleTypeTreeDatabase::Channel::gp_pPiPPiMPi0,
               makeBox(EPTrange,25.0) });
    db.insert({ParticleTypeTreeDatabase::Channel::gp_pPiPPiMPi0Pi0,
               makeBox(EPTrange,18.0) });

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
        // recover for below threshold
        if (energy < dataE.front()){
            return  dataXsec.front();
        }
        // don't crash for unknown energies, just use last known value
        if (energy > dataE.back()){
            LOG(WARNING) << "Interpolation for E = " << energy << " out of data-range, using last datapoint.";
            return  dataXsec.back();
        }

        auto xsec = ROOT::Math::Interpolator(dataE, dataXsec).Eval(energy);

        //fix for Interpolator smoothing into negative numbers
        if (xsec < 0 ) xsec = 0;

        return xsec;
    };
}

