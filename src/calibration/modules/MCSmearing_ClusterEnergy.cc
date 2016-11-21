#include "MCSmearing_ClusterEnergy.h"

#include "calibration/DataManager.h"

#include "tree/TCalibrationData.h"
#include "tree/TDetectorReadHit.h"
#include "tree/TCluster.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TH2.h"
#include "TH3.h"

#include <cstdint>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <vector>
#include <list>
#include <cmath>

#include "TRandom.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::std_ext;

MCSmearing_ClusterEnergy::MCSmearing_ClusterEnergy(Detector_t::Type_t detectorType,
                           std::shared_ptr<DataManager> calmgr
                           , std::vector<double> defaultEnergySigma) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_"
        << "MCSmearing_Energy"
           ),
    DetectorType(detectorType),
    calibrationManager(calmgr),
    EnergySigma(defaultEnergySigma,"Energy Sigmas")
{}

MCSmearing_ClusterEnergy::~MCSmearing_ClusterEnergy()
{
}

void MCSmearing_ClusterEnergy::ApplyTo(clusters_t& clusters)
{
    const auto& entry = clusters.find(DetectorType);

    if(entry != clusters.end()) {

        for(auto& cluster : entry->second) {
            const auto sigma = EnergySigma.Get(cluster.CentralElement, cluster.Energy);
            cluster.Energy = gRandom->Gaus(cluster.Energy, sigma);
        }
    }
}

double MCSmearing_ClusterEnergy::CalibType::Get(unsigned channel, const double E) const {
    //TODO: implement. Get from interpolator
    return 0.0;
}

std::list<Updateable_traits::Loader_t> MCSmearing_ClusterEnergy::GetLoaders()
{

    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {
            const string file = formatter() << calibrationManager->GetCalibrationDataFolder()
                                            << "/" << "MCSmearing_ClusterEnergy.root";
//            WrapTFileInput f(file);

//            TH2* h;
//            f.GetObject("Energy", h);

            LOG(INFO) << file;

            //TODO: load TH2 from root file and feed to interpolator


        }
    };
}


















