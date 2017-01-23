#include "ClippedInterpolatorWrapper.h"
#include "TH2D.h"
#include "base/Array2D.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"

using namespace std;
using namespace ant;


struct OrnderedGrid2D {
    std::vector<double> x;
    std::vector<double> y;
    ant::Array2D        z;

    OrnderedGrid2D(const unsigned nx, const unsigned ny, const double dflt=0.0):
        x(nx),
        y(ny),
        z(nx, ny, dflt) {}

};


ostream& ant::ClippedInterpolatorWrapper::Print(ostream& stream) const
{
    stream << "Theta Range: " << xrange << " E Range: " << yrange;
    return stream;
}

void ant::ClippedInterpolatorWrapper::setInterpolator(ClippedInterpolatorWrapper::interpolator_ptr_t i) {
    interp = move(i);
    xrange = interp->getXRange();
    yrange = interp->getYRange();
}

double ant::ClippedInterpolatorWrapper::boundsCheck_t::clip(double v) const
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

std::unique_ptr<const Interpolator2D> ClippedInterpolatorWrapper::makeInterpolator(TH2D* hist) {

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

    return std_ext::make_unique<Interpolator2D>(grid.x, grid.y, grid.z.Data());
}

ostream&ant::ClippedInterpolatorWrapper::boundsCheck_t::Print(ostream& stream) const
{
    stream << range << "-> [" << underflow << "|" << unclipped << "|" << overflow << "]";
    return stream;
}

double ant::ClippedInterpolatorWrapper::GetPoint(double x, double y) const
{
    x = xrange.clip(x);
    y = yrange.clip(y);
    return interp->GetPoint(x,y);
}

ant::ClippedInterpolatorWrapper::~ClippedInterpolatorWrapper()
{}

ant::ClippedInterpolatorWrapper::ClippedInterpolatorWrapper(): interp(nullptr), xrange({0,0}), yrange({0,0})
{}

ant::ClippedInterpolatorWrapper::ClippedInterpolatorWrapper(ClippedInterpolatorWrapper::interpolator_ptr_t i):
    interp(move(i)), xrange(interp->getXRange()), yrange(interp->getYRange())
{}
