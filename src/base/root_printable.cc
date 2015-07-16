#include "root_printable.h"
#include "TVector3.h"

std::ostream &operator<<(std::ostream &stream, const TVector3 &v)
{
    stream << "(" << v.X() << "," << v.Y() << "," << v.Z();
    return stream;
}
