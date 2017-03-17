#include "ProductionDataBase.h"

#include "base/interval.h"
#include "base/Logger.h"


#include "detail/gp_pg.h"
#include "detail/gp_pPi0.h"
#include "detail/gp_pPi0Pi0.h"
#include "detail/gp_pEta.h"
#include "detail/gp_pEtaPi0.h"
#include "detail/gp_pEtaPrime.h"
#include "detail/gp_pOmega.h"
#include "detail/gp_pRho.h"
#include "detail/gp_SigmaPlusK0S.h"

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
         gp_SigmaPlusK0S,
         gp_nPiP}
         )
    {
        db.insert(ch);
    }

    // Estimations for EPT-Range only !!!!
    IntervalD EPTrange({1410,1600});
    //sergey-estimates:
    db.insert({ParticleTypeTreeDatabase::Channel::gp_p3Pi0,
               makeBox(EPTrange,1.2) });

    // following 2 are first datapoint taken from:
    // @article{STRUCZINSKI197645,
    //    title = "Study of photoproduction on hydrogen in a streamer chamber with tagged photons for 1.6 GeV < Eγ < 6.3 GeV Topological and reaction cross sections",
    //    journal = "Nuclear Physics B",
    //    volume = "108",
    //    number = "1",
    //    pages = "45 - 74",
    //    year = "1976",
    //    note = "",
    //    issn = "0550-3213",
    //    doi = "http://dx.doi.org/10.1016/0550-3213(76)90123-1",
    //    url = "http://www.sciencedirect.com/science/article/pii/0550321376901231",
    //    author = "W. Struczinski and P. Dittmann and V. Eckardt and P. Joos and A. Ladage and H. Meyer and D. Notz and G. Hentschel and J. Knobloch and E. Rabe and H. Taureg and M. Grimm and I. Derado and P. Schacht and R. Meinke",
    //    }
    db.insert({ParticleTypeTreeDatabase::Channel::gp_pPiPPiMPi0,  // note: subtract already included channels
               makeBox(EPTrange,25.0 - 6.44) });                  //       taking BR of final states into account
    db.insert({ParticleTypeTreeDatabase::Channel::gp_pPiPPiMPi0Pi0,
               makeBox(EPTrange,18.0 - 0.9) });

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

