#pragma once

#include <vector>

namespace ant {
namespace calibration {

struct Helpers {
    /**
     * @brief FitXY fits the given data points (x,y,sigma_y) to a straight line f(x)=a+b*x
     * @param x_ input x positions
     * @param y_ input y positions
     * @param sigy standard error of y data points
     * @param a y intercept
     * @param b slope
     */
    static void FitXY(const std::vector<double>& x_,
                      const std::vector<double>& y_,
                      const std::vector<double>& sigy,
                      double& a, double& b) {
        double S   = 0;
        double Sx  = 0;
        double Sy  = 0;
        double Sxx = 0;
        double Sxy = 0;
        for(size_t i=0;i<y_.size();i++) {
            const double s = sigy[i];
            const double s2 = s*s;
            const double x  = x_[i];
            const double y  = y_[i];
            S   += 1/s2;
            Sx  += x/s2;
            Sy  += y/s2;
            Sxx += x*x/s2;
            Sxy += x*y/s2;
        }
        const double D = S*Sxx-Sx*Sx;
        a = (Sxx*Sy - Sx*Sxy)/D;
        b = (S*Sxy-Sx*Sy)/D;
    }
};

}
}