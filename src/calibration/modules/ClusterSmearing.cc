#include "ClusterSmearing.h"

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

ClusterSmearing::ClusterSmearing(std::shared_ptr<ClusterDetector_t> det,
                           std::shared_ptr<DataManager> calmgr
                           ) :
    Calibration::BaseModule(
        std_ext::formatter()
        << Detector_t::ToString(det->Type)
        << "_"
        << "ClusterSmearing"
           ),
    DetectorType(det->Type),
    nelements(det->GetNChannels()),
    calibrationManager(calmgr),
    interpolator(nullptr)
{}

ClusterSmearing::~ClusterSmearing()
{
}


struct ClusterSmearing::SigmaInterpolator {

    // TODO: insert interpolator here
    double GetSigma(const unsigned element, const double E) const {
        const auto bin = hist->FindBin(E, element);
        return hist->GetBinContent(bin);
    }

    SigmaInterpolator(TH2* h): hist(h) {}

    TH2* hist;
};

void ClusterSmearing::ApplyTo(clusters_t& clusters)
{
    // only run if smearing data present. (sould be nullptr for "data" -> skip)
    if(interpolator) {

        const auto& entry = clusters.find(DetectorType);

        if(entry != clusters.end()) {

            for(auto& cluster : entry->second) {
                const auto sigma = interpolator->GetSigma(cluster.CentralElement, cluster.Energy);
                cluster.Energy = gRandom->Gaus(cluster.Energy, sigma);
            }
        }
    }
}


std::list<Updateable_traits::Loader_t> ClusterSmearing::GetLoaders()
{

    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {

            auto obj = calibrationManager->GetTObject(GetName(), "energy", currPoint, nextChangePoint);

            TH2* hist = dynamic_cast<TH2*>(obj);

            if(hist) {
                if(unsigned(hist->GetNbinsY()) == nelements) {
                    VLOG(3) << "Smearing Histogram found";
                    this->interpolator = std_ext::make_unique<SigmaInterpolator>(hist);
                } else {
                    LOG(INFO) << "Smearing histogram found, but wrong number of elements! Not using this one. Deactivating Smearing.";
                    this->interpolator = nullptr;
                }
            } else {
                VLOG(3) << "No Smearing histogram found! Deactivating Smearing.";
                this->interpolator = nullptr;
            }

        }
    };
}


















