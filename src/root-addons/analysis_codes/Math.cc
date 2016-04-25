#include "Math.h"
#include "base/math_functions/AsymGaus.h"


TF1*ant::Math::AsymGaus()
{
    return math::AsymGaus::GetTF1();
}
