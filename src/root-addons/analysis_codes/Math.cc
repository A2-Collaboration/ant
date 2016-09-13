#include "Math.h"
#include "base/math_functions/AsymGaus.h"
#include "base/math_functions/CrystalBall.h"


TF1*ant::Math::AsymGaus()
{
    return math::AsymGaus::GetTF1();
}

TF1 *ant::Math::CrystalBall()
{
    return math::CrystalBall::GetTF1();
}
