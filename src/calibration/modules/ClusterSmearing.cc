#include "ClusterSmearing.h"

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
                           const std::string &Name,
                           std::shared_ptr<DataManager> calmgr
                           ) :
    Calibration::BaseModule(
        std_ext::formatter()
        << Detector_t::ToString(det->Type)
        << "_"
        << Name
           ),
    DetectorType(det->Type),
    calibrationManager(calmgr),
    interpolator(nullptr)
{}

ClusterCorrection::~ClusterCorrection()
{
}


struct ClusterCorrection::Interpolator {

    // TODO: insert interpolator here
    double GetSigma(const double E, const double theta) const {
        return interp->GetPoint(E, cos(theta));
    }

    Interpolator(TH2D* h): hist(h) { CleanupHistogram(hist); CreateInterpolators(hist); }

    TH2D* hist;

    std::unique_ptr<ClippedInterpolatorWrapper> interp;

    void CreateInterpolators(TH2D* hist) {
        interp = std_ext::make_unique<ClippedInterpolatorWrapper>(
                     ClippedInterpolatorWrapper::makeInterpolator(hist)
                     );
    }

    static void CleanupHistogram(TH2* hist) {

        auto check = [] (const double x) {
            return isfinite(x) && x >= 0.0;
        };

        for(int y = 1; y<=hist->GetNbinsY(); ++y) {
            for(int x = 1; x<=hist->GetNbinsX(); ++x) {
                if(!check(hist->GetBinContent(x,y))) {
                    for(int dx=1; dx<=hist->GetNbinsX();++dx) {
                        if(x-dx >= 1 && check(hist->GetBinContent(x-dx,y))) {
                            hist->SetBinContent(x,y, hist->GetBinContent(x-dx,y));
                            break;
                        }
                        if(x+dx <= hist->GetNbinsX() && check(hist->GetBinContent(x+dx,y))) {
                            hist->SetBinContent(x,y, hist->GetBinContent(x+dx,y));
                            break;
                        }
                    }
                }
            }
        }
    }
};

void ClusterCorrection::ApplyTo(clusters_t& clusters)
{

    if(interpolator) {

        const auto& entry = clusters.find(DetectorType);

        if(entry != clusters.end()) {

            for(auto& cluster : entry->second) {
                ApplyTo(cluster);
            }
        }
    }
}


std::list<Updateable_traits::Loader_t> ClusterCorrection::GetLoaders()
{

    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {

            TCalibrationData cdata;

            if(!calibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint)){
                LOG(WARNING) << "No data found for " << GetName();
                this->interpolator = nullptr;
                return;
            }

            auto hist = detail::TH2Storage::Decode(cdata);

            this->interpolator = std_ext::make_unique<Interpolator>(hist);

        }
    };
}

void ClusterSmearing::ApplyTo(TCluster &cluster)
{
    const auto sigma  = interpolator ?interpolator->GetSigma(cluster.Energy, cluster.Position.Theta()) : 0.0;
    cluster.Energy    = gRandom->Gaus(cluster.Energy, sigma);
}

void ClusterScaling::ApplyTo(TCluster &cluster)
{
    const auto factor  = interpolator ? interpolator->GetSigma(cluster.Energy, cluster.Position.Theta()) : 1.0;
    cluster.Energy    *= factor;
}
