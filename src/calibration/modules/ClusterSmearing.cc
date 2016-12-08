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
#include "base/Interpolator.h"
#include <ostream>
#include "base/Array2D.h"

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

//---------------- copy & paste from interpolated pulls ----------------------------

struct OrnderedGrid2D {
    vector<double> x;
    vector<double> y;
    Array2D        z;

    OrnderedGrid2D(const unsigned nx, const unsigned ny, const double dflt=0.0):
        x(nx),
        y(ny),
        z(nx, ny, dflt) {}

};


struct ClippedInterpolatorWrapper : ant::printable_traits {
    using interpolator_ptr_t = std::unique_ptr<const Interpolator2D>;

    interpolator_ptr_t interp;

    struct boundsCheck_t : ant::printable_traits {
        interval<double> range;
        mutable unsigned underflow = 0;
        mutable unsigned unclipped = 0;
        mutable unsigned overflow  = 0;

        double clip(double v) const;

        boundsCheck_t(const interval<double> r): range(r) {}

        std::ostream& Print(std::ostream& stream) const override;
    };

    boundsCheck_t xrange;
    boundsCheck_t yrange;

    ClippedInterpolatorWrapper(interpolator_ptr_t i);
    ClippedInterpolatorWrapper();
    ~ClippedInterpolatorWrapper();
    double GetPoint(double x, double y) const;

    void setInterpolator(interpolator_ptr_t i);

    std::ostream& Print(std::ostream& stream) const override;

};

ClippedInterpolatorWrapper::ClippedInterpolatorWrapper(ClippedInterpolatorWrapper::interpolator_ptr_t i):
    interp(move(i)), xrange(interp->getXRange()), yrange(interp->getYRange())
{}

ClippedInterpolatorWrapper::ClippedInterpolatorWrapper(): interp(nullptr), xrange({0,0}), yrange({0,0})
{}

ClippedInterpolatorWrapper::~ClippedInterpolatorWrapper()
{}

double ClippedInterpolatorWrapper::GetPoint(double x, double y) const
{
    x = xrange.clip(x);
    y = yrange.clip(y);
    return interp->GetPoint(x,y);
}

void ClippedInterpolatorWrapper::setInterpolator(ClippedInterpolatorWrapper::interpolator_ptr_t i) {
    interp = move(i);
    xrange = interp->getXRange();
    yrange = interp->getYRange();
}

ostream& ClippedInterpolatorWrapper::Print(ostream& stream) const
{
    stream << "Theta Range: " << xrange << " E Range: " << yrange;
    return stream;
}


double ClippedInterpolatorWrapper::boundsCheck_t::clip(double v) const
{
    if(v < range.Start()) {
        underflow++;
        return range.Start();
    }

    if(v > range.Stop()) {
        overflow++;
        return range.Stop();
    }

    unclipped++;

    return v;
}

ostream& ClippedInterpolatorWrapper::boundsCheck_t::Print(ostream& stream) const
{
    stream << range << "-> [" << underflow << "|" << unclipped << "|" << overflow << "]";
    return stream;
}

std::unique_ptr<const Interpolator2D> makeInterpolator(TH2D* hist) {

    const unsigned nx = unsigned(hist->GetNbinsX());
    const unsigned pad_x = nx > 3 ? 1 : 2;

    const unsigned ny = unsigned(hist->GetNbinsY());
    const unsigned pad_y = ny > 3 ? 1 : 2;

    OrnderedGrid2D grid(nx+pad_x*2, ny+pad_y*2, std_ext::NaN);

    // extend x bin positions
    {
        double avgBinWidth = .0;
        for(unsigned x=0; x<nx; ++x) {
            grid.x.at(x+pad_x) = hist->GetXaxis()->GetBinCenter(int(x+1));
            avgBinWidth += hist->GetXaxis()->GetBinWidth(int(x+1));
        }
        avgBinWidth /= nx;
        for(unsigned x=1; x<=pad_x;++x) {
            grid.x.at(pad_x-x)    = grid.x.at(pad_x-x+1)    - avgBinWidth;
            grid.x.at(pad_x+nx+x-1) = grid.x.at(pad_x+nx+x-2) + avgBinWidth;
        }
    }

    // eytend y bin positions
    {
        double avgBinWidth = .0;
        for(unsigned y=0; y<ny; ++y) {
            grid.y.at(y+pad_y) = hist->GetYaxis()->GetBinCenter(int(y+1));
            avgBinWidth += hist->GetYaxis()->GetBinWidth(int(y+1));
        }
        avgBinWidth /= ny;
        for(unsigned y=1; y<=pad_y;++y) {
            grid.y.at(pad_y-y)    = grid.y.at(pad_y-y+1)    - avgBinWidth;
            grid.y.at(pad_y+ny+y-1) = grid.y.at(pad_y+ny+y-2) + avgBinWidth;
        }
    }

    // copy data to the "middle". leave some borders around
    grid.z.CopyRect(Array2D_TH2D(hist), pad_x, pad_y);



    // fill borders:

    // top and botton rows
    for(unsigned y=1; y<=pad_y; ++y ) {
        for(unsigned x=0; x<nx; ++x) {
            grid.z.at(x+pad_x,pad_y-y)       = grid.z.at(x+pad_x,pad_y-y+1);
            grid.z.at(x+pad_x,pad_y+ny+y-1)    = grid.z.at(x+pad_x,pad_y+ny+y-2);
        }
    }

    // first and last column
    for(unsigned x=1; x<=pad_x; ++x ) {
        for(unsigned y=0; y<ny; ++y) {
            grid.z.at(pad_x-x,   y+pad_y)      = grid.z.at(pad_x-x+1,   y+pad_y);
            grid.z.at(pad_x+nx+x-1,y+pad_y)    = grid.z.at(pad_x+nx+x-2,y+pad_y);
        }
    }

    // flood fill averages the rest
    FloodFillAverages::fillNeighborAverages(grid.z);

    return std_ext::make_unique<Interpolator2D>(grid.x,grid.y, grid.z.Data());
}

//-----------------------------------------------------------------------------------
struct ClusterSmearing::SigmaInterpolator {

    // TODO: insert interpolator here
    double GetSigma(const double E, const double theta) const {
        return interp->GetPoint(E, cos(theta));
    }

    SigmaInterpolator(TH2D* h): hist(h) { CleanupHistogram(hist); CreateInterpolators(hist); }

    TH2D* hist;

    std::unique_ptr<const Interpolator2D> interp;

    void CreateInterpolators(TH2D* hist) {
        interp = makeInterpolator(hist);
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
                const auto sigma = interpolator->GetSigma(cluster.Energy, cluster.Position.Theta());
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

            TH2D* hist = dynamic_cast<TH2D*>(obj);

            if(hist) {
                VLOG(3) << "Smearing Histogram found";
                this->interpolator = std_ext::make_unique<SigmaInterpolator>(hist);
            } else {
                VLOG(3) << "No Smearing histogram found! Deactivating Smearing.";
                this->interpolator = nullptr;
            }

        }
    };
}


















