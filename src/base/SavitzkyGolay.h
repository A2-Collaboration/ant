#pragma once

#include "base/interval.h"

#include <memory>
#include <vector>

namespace ant {

struct SavitzkyGolay {

    // window = window_left + window_right + 1
    SavitzkyGolay(int window, int polynom_order);
    SavitzkyGolay(int window_left, int window_right, int polynom_order);

    std::vector<double> Smooth(const std::vector<double>& y) const;

    template<typename GetY, typename SetY>
    void Convolute(const GetY& getY, const SetY& setY,
                   const interval<int>& range) const {
        double convolution = 0.0;
        const int points = n_l + n_r + 1;
        for (int k = 0; k < points; k++) {
            auto i = k - n_l; // i runs from -n_l to n_r (inclusive), -n_l <= i <= n_r

            // do some wrap around to keep it in range
            if(i<range.Start())
                i = range.Start() + (range.Start() - i);
            else if(i>range.Stop())
                i = range.Stop()  - (i - range.Stop() );
            convolution += gsl_matrix_get(h, n_l, k) * getY(i);
        }
        setY(convolution); // implicitly assume i=0
    }

protected:
    const int n_l;
    const int n_r;
    const int m;

    template<typename T>
    using deleted_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

    struct gsl_matrix;
    using gsl_matrix_t = deleted_unique_ptr<gsl_matrix>;
    gsl_matrix_t h;

    static double gsl_matrix_get(const gsl_matrix_t& m, const size_t i, const size_t j);
};

}