#pragma once

#include "base/types.h"
#include "base/printable.h"
#include "base/std_ext/math.h"

#include "base/vec3.h"

#include <list>
#include <memory>

namespace ant {

struct TTarget : printable_traits {

    vec3 Vertex;

    TTarget(const vec3& vertex = {std_ext::NaN, std_ext::NaN, std_ext::NaN}):
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
