#include "Array2D.h"
#include <stdexcept>
#include "base/std_ext/string.h"
#include <base/interval.h>
#include "TH2D.h"

using namespace ant;
using namespace ant::std_ext;
using namespace std;

Array2D_TH2D::~Array2D_TH2D()
{}

void Array2D_TH2D::Set(const unsigned x, const unsigned y, const double v)
{
    if(x>=Width())
        throw std::out_of_range(formatter() << "X index out of range:" << x << " / " << Width());

    if(y>=Height())
        throw std::out_of_range(formatter() << "Y index out of range:" << y << " / " << Height());

    hist->SetBinContent(int(x+1), int(y+1), v);
}

double Array2D_TH2D::Get(const unsigned x, const unsigned y) const
{
    if(x>=Width())
        throw std::out_of_range(formatter() << "X index out of range:" << x << " / " << Width());

    if(y>=Height())
        throw std::out_of_range(formatter() << "Y index out of range:" << y << " / " << Height());

    return hist->GetBinContent(int(x+1),int(y+1));
}

unsigned Array2D_TH2D::Width() const
{
    return unsigned(hist->GetNbinsX());
}

unsigned Array2D_TH2D::Height() const
{
    return unsigned(hist->GetNbinsY());
}

size_t Array2D_TH2D::Size() const
{
    return Width() * Height();
}

Array2D::~Array2D() {}

double&Array2D::at(const unsigned x, const unsigned y) {

    if(x>=Width())
        throw std::out_of_range(formatter() << "X index out of range:" << x << " / " << Width());

    if(y>=Height())
        throw std::out_of_range(formatter() << "Y index out of range:" << y << " / " << Height());

    return data.at(Bin(x,y));
}

const double&Array2D::at(const unsigned x, const unsigned y) const {

    if(x>=Width())
        throw std::out_of_range(formatter() << "X index out of range:" << x << " / " << Width());

    if(y>=Height())
        throw std::out_of_range(formatter() << "Y index out of range:" << y << " / " << Height());

    return data.at(Bin(x,y));
}

void Array2D::Set(const unsigned x, const unsigned y, const double v) {

    if(x>=Width())
        throw std::out_of_range(formatter() << "X index out of range:" << x << " / " << Width());

    if(y>=Height())
        throw std::out_of_range(formatter() << "Y index out of range:" << y << " / " << Height());

    at(x,y) = v;
}

double Array2D::Get(const unsigned x, const unsigned y) const {

    if(x>=Width())
        throw std::out_of_range(formatter() << "X index out of range:" << x << " / " << Width());

    if(y>=Height())
        throw std::out_of_range(formatter() << "Y index out of range:" << y << " / " << Height());

    return at(x,y);
}

Array2DBase::~Array2DBase() {}



void Array2DBase::CopyRect(const ant::Array2DBase& src, const unsigned x, const unsigned y)
{
    auto dst_bounds = interval2D<unsigned>({0, this->Width()}, {0, this->Height()});
    auto target     = interval2D<unsigned>({x, x+src.Width()}, {y, y+src.Height()});
    target.Clip(dst_bounds);

    for(unsigned px=target.x.Start(); px<target.x.Stop(); ++px) {
        for(unsigned py=target.y.Start(); py<target.y.Stop(); ++py) {
            this->Set(px,py, src.Get(px-x, py-y));
        }
    }

}

void Array2DBase::CopyRect(const Array2DBase& src, const ant::interval2D<unsigned>& src_rect, const unsigned x, const unsigned y)
{
    const auto dst_bounds = interval2D<unsigned>({0, this->Width()}, {0, this->Height()});
    const auto src_bounds = interval2D<unsigned>({0, src.Width()}, {0, src.Height()});
    interval2D<unsigned> src_rect_ = src_rect;
    src_rect_.Clip(src_bounds);

    auto target = interval2D<unsigned>(src_rect_.x+x, src_rect_.y+y);
    target.Clip(dst_bounds);

    for(unsigned px=target.x.Start(); px<target.x.Stop(); ++px) {
        for(unsigned py=target.y.Start(); py<target.y.Stop(); ++py) {
            this->Set(px,py, src.Get(px-x, py-y));
        }
    }
}

double FloodFillAverages::getNeighborAverage(const Array2DBase& hist, const unsigned x, const unsigned y) {

    double sum = {};
    unsigned n = 0;

    for(const auto& d : neighbors4) {
        const int bx = int(x) + d.first;
        const int by = int(y) + d.second;

        if(IsBinValid(hist,bx,by)) {
            const auto v = hist.Get(bx,by);
            if(std::isfinite(v)) {
                sum += v;
                ++n;
            }
        }
    }

    return n>0 ? sum / n : 0.0;
}

unsigned FloodFillAverages::getNeighborCount(const Array2DBase& hist, const unsigned x, const unsigned y) {

    unsigned n = 0;

    for(const auto& d : neighbors4) {
        const int bx = x + d.first;
        const int by = y + d.second;

        const auto valid = IsBinValid(hist,bx,by);
        if( valid && isfinite(hist.Get(unsigned(bx),unsigned(by)))) {
            ++n;
        }
    }

    return n;
}

void FloodFillAverages::fillNeighborAverages(Array2DBase& hist) {


    unsigned neighbors=0;

    do {
        neighbors=0;
        unsigned p_x =0;
        unsigned p_y =0;

        for(unsigned x=0; x<hist.Width(); ++x) {
            for(unsigned y=0; y<hist.Height(); ++y) {

                if(std::isnan(hist.Get(x,y))) {
                    const auto n = getNeighborCount(hist, x, y);
                    if(n>neighbors) {
                        neighbors = n;
                        p_x = x;
                        p_y = y;
                    }
                }
            }
        }

        // if updatable bin found
        if(neighbors > 0) {
            const auto a = getNeighborAverage(hist,p_x,p_y);
            hist.Set(p_x,p_y, a);
        }
    } while(neighbors != 0);
}

bool FloodFillAverages::IsBinValid(const Array2DBase& hist, int x, int y)
{
    return (x>=0) && (y>=0)
            && (unsigned(x)<hist.Width())
            && (unsigned(y)<hist.Height())
            && (isfinite(hist.Get(unsigned(x),unsigned(y))));
}

// up, down, right, left
const std::vector<std::pair<int,int>> FloodFillAverages::neighbors4 = {{+1,0},{-1,0},{0,+1},{0,-1}};
