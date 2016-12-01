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
#include "Math/Interpolator.h"

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
        return elements.at(element)->Eval(E);
    }

    SigmaInterpolator(TH2* h): hist(h) { CleanupHistogram(hist); CreateInterpolators(hist); }

    TH2* hist;

    std::vector<std::unique_ptr<ROOT::Math::Interpolator>> elements;

    void CreateInterpolators(const TH2* hist) {
        const auto bins = hist->GetNbinsX();
        const auto bw = hist->GetXaxis()->GetBinWidth(1);

        elements.reserve(unsigned(hist->GetNbinsY()));

        for(int y=1; y<=hist->GetNbinsY(); ++y) {

            vector<double> px; px.reserve(unsigned(bins+4));
            vector<double> py; py.reserve(unsigned(bins+4));

            px.push_back(hist->GetXaxis()->GetBinCenter(1) - 2*bw);
            py.push_back(hist->GetBinContent(1,y));

            px.push_back(hist->GetXaxis()->GetBinCenter(1) - 1*bw);
            py.push_back(hist->GetBinContent(1,y));

            for(int x=1; x<=hist->GetNbinsX(); ++x) {
                px.push_back(hist->GetXaxis()->GetBinCenter(x));
                py.push_back(hist->GetBinContent(x,y));
            }

            px.push_back(hist->GetXaxis()->GetBinCenter(bins) + 1*bw);
            py.push_back(hist->GetBinContent(bins,y));
            px.push_back(hist->GetXaxis()->GetBinCenter(bins) + 2*bw);
            py.push_back(hist->GetBinContent(bins,y));

            elements.emplace_back(std_ext::make_unique<ROOT::Math::Interpolator>(px,py));
        }
    }

    static void CleanupHistogram(TH2* hist) {
        for(int y = 1; y<=hist->GetNbinsY(); ++y) {
            for(int x = 1; x<=hist->GetNbinsX(); ++x) {
                if(hist->GetBinContent(x,y) < 0.0) {
                    for(int dx=1; dx<=hist->GetNbinsX();++dx) {
                        if(x-dx >= 1 && hist->GetBinContent(x-dx,y) >= 0.0) {
                            hist->SetBinContent(x,y, hist->GetBinContent(x-dx,y));
                            break;
                        }
                        if(x+dx <= hist->GetNbinsX() && hist->GetBinContent(x+dx,y) >= 0.0) {
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

            auto obj = calibrationManager->GetTObject(GetName(), "energy_smearing", currPoint, nextChangePoint);

            TH2* hist = dynamic_cast<TH2*>(obj);

            if(hist) {
                if(unsigned(hist->GetNbinsY()) == nelements) {
                    VLOG(3) << "Smearing Histogram found";
                    this->interpolator = std_ext::make_unique<SigmaInterpolator>(hist);
                } else {
                    LOG(WARNING) << "Smearing histogram found, but wrong number of elements! Not using this one. Deactivating Smearing.";
                    this->interpolator = nullptr;
                }
            } else {
                VLOG(3) << "No Smearing histogram found! Deactivating Smearing.";
                this->interpolator = nullptr;
            }

        }
    };
}


















