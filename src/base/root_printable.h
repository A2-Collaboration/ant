#ifndef ANT_ROOT_PRINTABLE_H
#define ANT_ROOT_PRINTABLE_H

#include <ostream>
#include "TVector3.h"

std::ostream &operator<<(std::ostream& stream, const TVector3& v);

#endif
