#include "ClusterCorrection.h"

#include "calibration/DataManager.h"
#include "tree/TCalibrationData.h"
#include "tree/TDetectorReadHit.h"
#include "tree/TCluster.h"
#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/ClippedInterpolatorWrapper.h"
#include "detail/TH2Storage.h"

#include "TH2.h"

#include <list>
#include <cmath>

#include "TRandom.h"



using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::std_ext;

ClusterCorrection::ClusterCorrection(std::shared_ptr<ClusterDetector_t> det,
                           const std::string &Name, const Filter_t Filter,
                           std::shared_ptr<DataManager> calmgr
                           ) :
    Calibration::BaseModule(
        std_ext::formatter()
        << Detector_t::ToString(det->Type)
        << "_"
        << Name
           ),
    DetectorType(det->Type),
    filter(Filter),
    calibrationManager(calmgr),
    interpolator(nullptr)
{}

ClusterCorrection::~ClusterCorrection()
{
}


void ClusterCorrection::ApplyTo(clusters_t& clusters)
{

    if(interpolator) {

        const auto& entry = clusters.find(DetectorType);

        if(entry != clusters.end()) {

            for(auto& cluster : entry->second) {

                ApplyTo(cluster);

                if(cluster.Energy < 0.0)
                    cluster.Energy = 0.0;
            }
        }
    }
}


std::list<Updateable_traits::Loader_t> ClusterCorrection::GetLoaders()
{

    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {

            const bool isMC    = currPoint.isSet(TID::Flags_t::MC);

            TCalibrationData cdata;

            const bool loadOK = calibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint);

            if((isMC && filter == Filter_t::Data)
                  || (!isMC  && filter == Filter_t::MC)) {
                this->interpolator = nullptr;
                return;
            }

            if(!loadOK) {
                LOG(WARNING) << "No data found for " << GetName();
                this->interpolator = nullptr;
                return;
            }

            auto hist = detail::TH2Storage::Decode(cdata);

            this->interpolator = std_ext::make_unique<ClippedInterpolatorWrapper>(
                                     ClippedInterpolatorWrapper::makeInterpolator(hist));

            delete hist;
        }
    };
}

void ClusterSmearing::ApplyTo(TCluster& cluster)
{
    const auto sigma  = interpolator->GetPoint(cluster.Energy, cos(cluster.Position.Theta()));
    cluster.Energy    = gRandom->Gaus(cluster.Energy, sigma);
}

void ClusterECorr::ApplyTo(TCluster& cluster)
{
    const auto factor  = interpolator->GetPoint(cluster.Energy, cluster.Hits.size());
    cluster.Energy    *= factor;
}
