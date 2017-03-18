#pragma once

#include "base/types.h"
#include "base/std_ext/math.h"

#include "base/vec/vec3.h"

#include <list>
#include <memory>

namespace ant {

struct TTarget {

    vec3 Vertex;

    TTarget(const vec3& vertex = {std_ext::NaN, std_ext::NaN, std_ext::NaN}):
        Vertex(vertex)
    {}


    template<class Archive>
    void serialize(Archive& archive) {
        archive(Vertex);
    }


    friend std::ostream& operator<<(std::ostream& s, const TTarget& o) {
        s << "Target Vertex=" << o.Vertex;
        return s;
    }

};

}
