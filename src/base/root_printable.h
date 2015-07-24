#pragma once

#include <ostream>
#include "TVector3.h"

inline std::ostream &operator<<(std::ostream& stream, const TVector3& v) {
    stream << "(" << v.X() << "," << v.Y() << "," << v.Z() << ")";
    return stream;
}


