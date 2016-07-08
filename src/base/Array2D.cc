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

template <typename T>
struct interval2D {
    interval<T> x;
    interval<T> y;

    interval2D(const interval<T>& xi, const interval<T>& yi): x(xi), y(yi) {}

    void Clip(const interval2D<T>& other) {
        x = intersect(x, other.x);
        y = intersect(y, other.y);
    }
};

void CopyRect(const Array2DBase& src, Array2DBase& dst, const unsigned x, const unsigned y)
{
    auto dst_bounds = interval2D<unsigned>({0, dst.Width()}, {0, dst.Height()});
    auto target     = interval2D<unsigned>({x, x+src.Width()}, {y, y+src.Height()});
    target.Clip(dst_bounds);

}
