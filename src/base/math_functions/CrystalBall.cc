#include "CrystalBall.h"

#include "TF1.h"

#include <cmath>

using namespace ant;
using namespace ant::math;

/**
 * @brief crystalball_function
 * @param x
 * @param alpha
 * @param n
 * @param sigma
 * @param mean
 * @return
 * @note copied from ROOT6 (https://root.cern.ch/doc/master/PdfFuncMathCore_8cxx_source.html), with little improvements
 */
double CrystalBall::Eval(const double x, const double alpha, const double n, const double sigma, const double mean, const double hight) noexcept {

    // evaluate the crystal ball function
    if (sigma < 0.0)     return 0.0;

    const double z = (x - mean) / sigma * ((alpha < 0.0) ? 1.0 : -1.0);

    const double abs_alpha = std::abs(alpha);

    if (z  > - abs_alpha)
        return std::exp(- 0.5 * z * z) * hight;

    else {
        const double nDivAlpha = n/abs_alpha;
        const double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
        const double B = nDivAlpha -abs_alpha;
        const double arg = nDivAlpha/(B-z);
        return AA * std::pow(arg,n) * hight;
    }
}

double CrystalBall::Eval_ROOT(const double *x, const double *p) {
    return Eval(x[0], p[0], p[1], p[2], p[3], p[4]);
}

TF1* CrystalBall::GetTF1() {
    auto f = new TF1("CrystalBall", Eval_ROOT, -10, 10, 5);
    f->SetParName(0, "#alpha");
    f->SetParName(1, "N");
    f->SetParName(2, "#sigma");
    f->SetParName(3, "Mean");
    f->SetParName(4, "Hight");
    return f;
}
