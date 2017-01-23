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
    calibrationManager(calmgr),
    smearing_interpolator(nullptr)
{}

ClusterSmearing::~ClusterSmearing()
{
}


struct ClusterSmearing::SigmaInterpolator {

    // TODO: insert interpolator here
    double GetSigma(const double E, const double theta) const {
        return interp->GetPoint(E, cos(theta));
    }

    SigmaInterpolator(TH2D* h): hist(h) { CleanupHistogram(hist); CreateInterpolators(hist); }

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

void ClusterSmearing::ApplyTo(clusters_t& clusters)
{

    // only run if smearing / scaling data present. (sould be nullptr for "data" -> skip)
    if(smearing_interpolator || smearing_interpolator) {

        const auto& entry = clusters.find(DetectorType);

        if(entry != clusters.end()) {

            for(auto& cluster : entry->second) {
                const auto factor = scaling_interpolator  ?  scaling_interpolator->GetSigma(cluster.Energy, cluster.Position.Theta()) : 1.0;
                const auto sigma  = smearing_interpolator ? smearing_interpolator->GetSigma(cluster.Energy, cluster.Position.Theta()) : 0.0;
                cluster.Energy = gRandom->Gaus(cluster.Energy * factor, sigma);
            }
        }
    }
}


std::list<Updateable_traits::Loader_t> ClusterSmearing::GetLoaders()
{

    return {
        [this] (const TID& currPoint, TID& nextChangePoint) {

            TCalibrationData cdata;

            if(!calibrationManager->GetData(GetName(), currPoint, cdata, nextChangePoint)){
                VLOG(3) << "No Cluster Smearings found";
                this->smearing_interpolator = nullptr;
                return;
            }

            auto hist = detail::TH2Storage::Decode(cdata);

            this->smearing_interpolator = std_ext::make_unique<SigmaInterpolator>(hist);

        },
        [this] (const TID& currPoint, TID& nextChangePoint) {

            TCalibrationData cdata;

            if(!calibrationManager->GetData(formatter() << GetName() << "_scaleing", currPoint, cdata, nextChangePoint)){
                VLOG(3) << "No Cluster Energy Scaling factor found";
                this->scaling_interpolator = nullptr;
                return;
            }

            auto hist = detail::TH2Storage::Decode(cdata);

            this->smearing_interpolator = std_ext::make_unique<SigmaInterpolator>(hist);

        }
    };
}
