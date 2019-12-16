#pragma once

#include <utility>

class TF1;

namespace ant {
namespace math {

struct SimpleCrystalBall {
    static double Eval(double x, double mean, double sigma, double k, double height) noexcept;

    static double Eval_ROOT(const double* x, const double* p);

    static TF1* GetTF1();

};

// alias for SimpleCrystalBall function
struct GaussExp {
    template <typename... Args>
    static auto GetTF1(Args&&... args) -> decltype(SimpleCrystalBall::GetTF1(std::forward<Args>(args)...))
    {
        return SimpleCrystalBall::GetTF1(std::forward<Args>(args)...);
    }
};

}
}
