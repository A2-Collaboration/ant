#include "Array2D.h"

#include "base/std_ext/string.h"
#include "base/interval.h"
#include "base/FloodFillAverages.h"
#include "base/std_ext/math.h"

#include "TH2D.h"

#include <stdexcept>

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

double& Array2D::at(const unsigned x, const unsigned y) {

    if(x>=Width())
        throw std::out_of_range(formatter() << "X index out of range:" << x << " / " << Width());

    if(y>=Height())
        throw std::out_of_range(formatter() << "Y index out of range:" << y << " / " << Height());

    return data.at(Bin(x,y));
}

const double& Array2D::at(const unsigned x, const unsigned y) const {

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

void Array2DBase::FloodFillAverages()
{
    auto getVal = [this] (int i) {
        return Get(i % Width(), i / Width());
    };
    auto setVal = [this] (int i, double v) {
        Set(i % Width(), i / Width(), v);
    };
    auto getNeighbours = [this] (int i) {
        const int x = i % Width();
        const int y = i / Width();
        vector<int> neighbours;
        for(auto& dir : vector<pair<int,int>>{{+1,0},{-1,0},{0,+1},{0,-1}}) {
            const auto bx = x + dir.first;
            const auto by = y + dir.second;
            if(bx<0 || bx >= int(Width()))
                continue;
            if(by<0 || by >= int(Height()))
                continue;
            neighbours.push_back(bx + by*Width());
        }
        return neighbours;
    };
    auto getValid = [getVal] (int i) {
        return isfinite(getVal(i));
    };

    floodFillAverages(Size(), getVal, setVal, getNeighbours, getValid);
}

void Array2DBase::RemoveOutliers(double IQR_factor_lo, double IQR_factor_hi)
{
    std_ext::IQR iqr;
    for(auto x=0u;x<Width();x++) {
        for(auto y=0u;y<Height();y++) {
            iqr.Add(Get(x,y));
        }
    }

    if(iqr.GetN()==0)
        return;

    const interval<double> valid_range(
                iqr.GetMedian() - IQR_factor_lo*iqr.GetIQR(),
                iqr.GetMedian() + IQR_factor_hi*iqr.GetIQR()
                );

    for(auto x=0u;x<Width();x++) {
        for(auto y=0u;y<Height();y++) {
            if(!valid_range.Contains(Get(x,y)))
                Set(x, y, std_ext::NaN);
        }
    }
}
