#include "SimpleCrystalBall.h"

#include "TF1.h"

#include <cmath>

using namespace ant;
using namespace ant::math;

/**
 * @brief Simple Crystal Ball Function
 * @param x
 * @param mean
 * @param sigma
 * @param k (decay constant)
 * @param height
 * @note
 * @return Function implemented based on https://arxiv.org/pdf/1603.08591.pdf
 */
double SimpleCrystalBall::Eval(double x, double mean, double sigma, double k, double height) noexcept {

    // evaluate the simple crystal ball function
    if (sigma < 0.0) return 0.0;

    const double z = (x - mean) / sigma;

    const double abs_k = std::abs(k);

    if (z >= -abs_k)
        return std::exp(-0.5*z*z) * height;
    else
        return std::exp(0.5*k*k + abs_k*z) * height;
}

double SimpleCrystalBall::Eval_ROOT(const double *x, const double *p) {
    return Eval(x[0], p[0], p[1], p[2], p[3]);
}

TF1* SimpleCrystalBall::GetTF1() {
    auto f = new TF1("SimpleCrystalBall", Eval_ROOT, -10, 10, 5);
    f->SetParName(0, "Mean");
    f->SetParName(1, "#sigma");
    f->SetParName(2, "Decay constant k");
    f->SetParName(3, "Height");
    return f;
}
