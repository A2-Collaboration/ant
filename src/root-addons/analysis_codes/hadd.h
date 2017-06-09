#pragma once

#include "TDirectory.h"
#include <memory>
#include <vector>

namespace ant {

struct hadd {

    template<typename T>
    using unique_ptrs_t = std::vector<std::unique_ptr<T>>;
    using sources_t = unique_ptrs_t<const TDirectory>;

    static void MergeRecursive(TDirectory& target, const sources_t& sources, unsigned& nPaths);

};

}
