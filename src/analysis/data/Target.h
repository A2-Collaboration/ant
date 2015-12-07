#pragma once

#include "base/types.h"
#include "base/printable.h"
#include "tree/TDataRecord.h"

#include "TVector3.h"
#include <list>
#include <memory>

namespace ant {
namespace analysis {
namespace data {

struct Target_t : printable_traits {

    Target_t( const TVector3& vertex = TVector3()):
        Vertex(vertex)
    {}

    TVector3 Vertex;

    std::ostream& Print(std::ostream& stream) const;
};

}
}
}
