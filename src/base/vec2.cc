#include "vec2.h"

#include "TVector2.h"

using namespace ant;


vec2::vec2(const TVector2& v) noexcept :
    x(v.X()),y(v.Y())
{}

vec2& vec2::operator=(const TVector2& v) noexcept {
    x = v.X();
    y = v.Y();
    return *this;
}

vec2::operator TVector2() const noexcept {
    return TVector2(x,y);
}
