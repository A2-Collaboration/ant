#include "MCSmearing_ClusterEnergy.h"

#include "calibration/DataManager.h"
#include "expconfig/detectors/CB.h"

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

MCSmearing_ClusterEnergy::MCSmearing_ClusterEnergy(std::shared_ptr<expconfig::detector::CB> cb,
                           std::shared_ptr<DataManager> calmgr
                           ) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(cb->Type)
        << "_"
        << "MCSmearing_Energy"
           ),
    DetectorType(cb->Type),
    calibrationManager(calmgr),
    EnergySigma({},"Energy Sigmas")
{}

MCSmearing_ClusterEnergy::~MCSmearing_ClusterEnergy()
{
}

void MCSmearing_ClusterEnergy::GetGUIs(std::list<std::unique_ptr<gui::CalibModule_traits> >&, const OptionsPtr) {}

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

            auto obj = calibrationManager->GetTObject(GetName(), "smearings", currPoint, nextChangePoint);

            TH2* hist = dynamic_cast<TH2*>(obj);

            if(hist) {
                LOG(INFO) << "Histogram found: ";
            }


            //TODO: load TH2 from root file and feed to interpolator


        }
    };
}


















