#include "SavitzkyGolay.h"

#include "base/std_ext/string.h"

extern "C" {
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
}

using namespace std;
using namespace ant;

// little trick to forward declare gsl_matrix, which is a POD
struct SavitzkyGolay::gsl_matrix : ::gsl_matrix {};

SavitzkyGolay::SavitzkyGolay(int window, int polynom_order) :
    // window = n_l + n_r + 1 = number of points
    // need to differentiate odd/even cases, prefer looking "back" (n_l>=n_r) if even case
    SavitzkyGolay((window-1)/2 + (window % 2 == 0), (window-1)/2, polynom_order)
{
}

SavitzkyGolay::SavitzkyGolay(int window_left, int window_right, int polynom_order) :
    n_l(window_left),
    n_r(window_right),
    m(polynom_order),
    h(MakeH(n_l,n_r,m))
{
}

SavitzkyGolay::gsl_unique_ptr<SavitzkyGolay::gsl_matrix> SavitzkyGolay::MakeH(int n_l, int n_r, int m)
{
    const auto points = n_l + n_r + 1;

    // define some unique_ptr alloc for matrix
    auto gsl_matrix_alloc = [] (size_t n1, size_t n2) {
        return gsl_unique_ptr<gsl_matrix>(static_cast<gsl_matrix*>(::gsl_matrix_alloc(n1, n2)), ::gsl_matrix_free);
    };

    // the code in the following will  eventually set h (the precomputed Savitzky Golay coefficients)
    auto h = gsl_matrix_alloc(points, points);

    // gsl_permutation is only used here, so provide just the exception-safe unique_ptr allocator
    auto gsl_permutation_alloc = [] (size_t n) {
        return gsl_unique_ptr<gsl_permutation>(static_cast<gsl_permutation*>(::gsl_permutation_alloc(n)), ::gsl_permutation_free);
    };

    // provide error handling, as we're using unique_ptr, the code is exception safe!
    auto catch_gsl_error = [] (int error) {
        if(error == GSL_SUCCESS)
            return;
        throw Exception(std_ext::formatter() << "Encountered GSL error: " << gsl_strerror(error));
    };

    // compute Vandermonde matrix
    auto vandermonde = gsl_matrix_alloc(points, m + 1);
    for (int i = 0; i < points; ++i){
        gsl_matrix_set(vandermonde, i, 0, 1.0);
        for (int j = 1; j <= m; ++j)
            gsl_matrix_set(vandermonde, i, j, gsl_matrix_get(vandermonde, i, j - 1) * i);
    }

    // compute V^TV
    auto vtv = gsl_matrix_alloc(m + 1, m + 1);
    catch_gsl_error(gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, vandermonde, vandermonde, 0.0, vtv));

    // compute (V^TV)^(-1) using LU decomposition
    auto p = gsl_permutation_alloc(m + 1);
    int signum;
    catch_gsl_error(gsl_linalg_LU_decomp(vtv, p, &signum));

    auto vtv_inv = gsl_matrix_alloc(m + 1, m + 1);
    catch_gsl_error(gsl_linalg_LU_invert(vtv, p, vtv_inv));

    // compute (V^TV)^(-1)V^T
    auto vtv_inv_vt = gsl_matrix_alloc(m + 1, points);
    catch_gsl_error(gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, vtv_inv, vandermonde, 0.0, vtv_inv_vt));

    // finally, compute H = V(V^TV)^(-1)V^T
    catch_gsl_error(gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, vandermonde, vtv_inv_vt, 0.0, h));

    return h;
}

vector<double> SavitzkyGolay::Smooth(const vector<double>& y) const
{
    const int points = n_l + n_r + 1;

    // for edge cases, mirror the input y as follows:
    // ... y[2] y[1] y[0] y[1] .... y[n-2] y[n-1] y[n] y[n+1] ...
    auto get_y = [&y] (const int i) {
        if(i<0)
            return y.at(-i);
        if(unsigned(i)>=y.size()) {
            const auto i_ = y.size()-(i-y.size())-2;
            return y.at(i_);
        }
        return y[i];
    };

    const int d_n = y.size();
    vector<double> result(d_n);
    for (int i = 0; i < d_n; i++){
        double convolution = 0.0;
        for (int k = 0; k < points; k++)
            convolution += gsl_matrix_get(h, n_l, k) * get_y(i - n_l + k);
        result[i] = convolution;
    }

    return result;
}

double SavitzkyGolay::gsl_matrix_get(const gsl_matrix* m, const size_t i, const size_t j)
{
    return ::gsl_matrix_get(m, i, j);
}

