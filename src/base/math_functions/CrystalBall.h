#pragma once

class TF1;

namespace ant {
namespace math {

struct CrystalBall {
    static double Eval(double x, double alpha, double n, double sigma, double mean, double hight) noexcept;

    static double Eval_ROOT(const double* x, const double* p);

    static TF1* GetTF1();

};

}
}
