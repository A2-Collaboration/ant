#pragma once

#include "base/types.h"
#include "base/printable.h"
#include "base/root_printable.h"
#include "base/std_ext/math.h"

#include "TVector3.h"

#include <list>
#include <memory>

namespace ant {

struct TTarget : printable_traits {

    TVector3 Vertex;

    TTarget(const TVector3& vertex = TVector3(std_ext::NaN, std_ext::NaN, std_ext::NaN)):
        Vertex(vertex)
    {}


    template<class Archive>
    void serialize(Archive& archive) {
        archive(Vertex);
    }


    std::ostream& Print(std::ostream& s) const {
        s << "Target Vertex=" << Vertex;
        return s;
    }

};

}
