#include "Interpolator.h"

extern "C" {
#include "detail/interp2d/interp2d_spline.h"
}

using namespace std;
using namespace ant;

struct Interpolator2D::interp2d : ::interp2d {};
struct Interpolator2D::gsl_interp_accel : ::gsl_interp_accel {};

const interp2d_type* getType(Interpolator2D::Type type) {
    switch(type) {
    case Interpolator2D::Type::Bilinear:
        return interp2d_bilinear;
    case Interpolator2D::Type::Bicubic:
        return interp2d_bicubic;
    }
    return nullptr;
}


Interpolator2D::Interpolator2D(const std::vector<double>& x,
                               const std::vector<double>& y,
                               const std::vector<double>& z,
                               Type type) :
    X(x), Y(y), Z(z),
    interp(static_cast<interp2d*>(
               interp2d_alloc(getType(type), X.size(), Y.size())
               ), interp2d_free),
    xa(static_cast<gsl_interp_accel*>(gsl_interp_accel_alloc()), gsl_interp_accel_free),
    ya(static_cast<gsl_interp_accel*>(gsl_interp_accel_alloc()), gsl_interp_accel_free)
{
    if(X.size()*Y.size() != Z.size())
        throw Exception("X*Y grid must match to Z values");
    interp2d_init(interp.get(), X.data(), Y.data(), Z.data(), X.size(), Y.size());
}

double Interpolator2D::GetPoint(double x, double y) const
{
    return interp2d_eval(interp.get(), X.data(), Y.data(), Z.data(), x, y,
                         xa.get(), ya.get());
}
