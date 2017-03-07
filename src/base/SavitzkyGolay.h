#pragma once

#include <memory>
#include <vector>

namespace ant {

struct SavitzkyGolay {

    SavitzkyGolay(int n_points_left, int n_points_right, int polynom_order);

    std::vector<double> Smooth(const std::vector<double>& y);

protected:
    const int n_l;
    const int n_r;
    const int m;

    template<typename T>
    using deleted_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

    struct gsl_matrix;
    using gsl_matrix_t = deleted_unique_ptr<gsl_matrix>;
    gsl_matrix_t h;
};

}