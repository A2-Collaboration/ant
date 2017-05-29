#pragma once

#include <vector>
#include <stdexcept>
#include <memory>
#include <functional>
#include "base/interval.h"

namespace ant {

class Interpolator2D {
public:
    enum class Type {
        Bilinear, Bicubic
    };

    Interpolator2D(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& z,
                   Type type = Type::Bicubic);

    double GetPoint(double x, double y) const;

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

    interval<double> getXRange() const;
    interval<double> getYRange() const;


private:

    const std::vector<double> X;
    const std::vector<double> Y;
    const std::vector<double> Z;

    template<typename T>
    using deleted_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

    struct interp2d;
    deleted_unique_ptr<interp2d> interp;

    struct gsl_interp_accel;
    deleted_unique_ptr<gsl_interp_accel> xa;
    deleted_unique_ptr<gsl_interp_accel> ya;
};

}
