#include "Interpolator.h"

extern "C" {
#include "detail/interp2d/interp2d_spline.h"
}

using namespace std;
using namespace ant;


Interpolator2D::Interpolator2D(const std::vector<double>& x,
                               const std::vector<double>& y,
                               const std::vector<double>& z) :
    X(x), Y(y), Z(z)
{
    if(X.size()*Y.size() != Z.size())
        throw Exception("X*Y grid must match to Z values");

    interp2d_spline* spline = interp2d_spline_alloc(interp2d_bicubic, X.size(), Y.size());

    (void)spline;
}
