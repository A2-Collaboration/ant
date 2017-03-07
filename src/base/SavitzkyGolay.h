#pragma once

#include <memory>
#include <vector>

namespace ant {

struct SavitzkyGolay {

    SavitzkyGolay(int n_l, int n_r, int m);

    std::vector<double> Smooth(const std::vector<double>& y);

protected:
    int d_sav_gol_points = 2;
    int d_smooth_points = 2;
    int d_polynom_order = 2;

    template<typename T>
    using deleted_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

    struct gsl_matrix;
    using gsl_matrix_t = deleted_unique_ptr<gsl_matrix>;
    gsl_matrix_t h;
};

}