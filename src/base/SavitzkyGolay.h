#pragma once

#include "base/interval.h"

#include <memory>
#include <vector>
#include <functional>

namespace ant {

struct SavitzkyGolay {

    // window = window_left + window_right + 1
    SavitzkyGolay(int window, int polynom_order);
    SavitzkyGolay(int window_left, int window_right, int polynom_order);

    std::vector<double> Smooth(const std::vector<double>& y) const;

    template<typename GetY, typename SetY>
    void Convolute(const GetY& getY, const SetY& setY,
                   const interval<int>& range) const
    {
        double convolution = 0.0;
        const auto points = n_l + n_r + 1;
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

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

protected:
    const int n_l;
    const int n_r;
    const int m;

    template<typename T, typename Deleter = std::function<void(T*)>>
    struct gsl_unique_ptr : std::unique_ptr<T, Deleter>
    {
        using std::unique_ptr<T, Deleter>::unique_ptr;
        // this implicit conversion makes the type transparently usable with gsl_* methods
        operator T* () const { return this->get(); }
    };

    struct gsl_matrix;
    const gsl_unique_ptr<gsl_matrix> h;
    static gsl_unique_ptr<gsl_matrix> MakeH(int n_l, int n_r, int m);

    // define wrapper for above templated Convolute method
    static double gsl_matrix_get(const gsl_matrix* m, const size_t i, const size_t j);
};

}