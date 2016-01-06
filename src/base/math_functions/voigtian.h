#pragma once

class TF1;

namespace ant {
namespace math {

/**
 * @brief Voigtian profile
 * @see https://en.wikipedia.org/wiki/Voigt_profile
 */
struct voigtian {

    /**
     * @brief Evaluate
     * @param x
     * @param A coefficient . Not amplitude -> not normalized
     * @param mu position of the maximum
     * @param sigma
     * @param gamma
     * @return
     */
    static double Eval(const double x, const double A, const double mu, const double sigma, const double gamma);

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
