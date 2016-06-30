#pragma once

#include <vector>
#include <stdexcept>

namespace ant {

class Interpolator2D {
public:
    Interpolator2D(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& z);

    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

protected:
    const std::vector<double> X;
    const std::vector<double> Y;
    const std::vector<double> Z;

};

}