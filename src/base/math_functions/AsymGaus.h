#pragma once

class TF1;

namespace ant {
namespace math {

/**
 * @brief Asymmetric Gaus
 * Left and right side of mu have different sigmas
 */
struct AsymGaus {

    /**
     * @brief Evaluate
     * @param x
     * @param A Amplitude
     * @param mu position of the maximum
     * @param sigma_l
     * @param sigma_r
     * @return
     */
    static double Eval(const double x, const double A, const double mu, const double sigma_l, const double sigma_r);

    /**
     * @brief Evaluate (ROOT style)
     * @param x ptr to x variables (only element 0 used)
     * @param p ptr to parameters (0..3 used)
     * @return
     */
    static double Eval_ROOT(const double* x, const double* p);

    static TF1* GetTF1();
};
}
}
