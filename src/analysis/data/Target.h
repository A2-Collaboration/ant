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

class Target_t : public printable_traits {
protected:

    TVector3 vertex;

public:

    Target_t( const TVector3& vertex_ = TVector3()):
        vertex(vertex_)
    {}

    TVector3  Vertex() const { return vertex; }
    TVector3& Vertex()       { return vertex; }

    std::ostream& Print(std::ostream& stream) const;
};

}
}
}
