#pragma once

#include <vector>
#include <algorithm>
#include <numeric>

namespace ant {
namespace std_ext {


/**
 * @brief Takes a vector and compares its contents to obtain sorted indices
 * @param vec vector containing comparable objects
 * @return vector<size_t> of sorted indices
 */
template<typename T>
std::vector<size_t> get_sorted_indices(std::vector<T> vec)
{
    std::vector<size_t> v(vec.size());
    std::iota(v.begin(), v.end(), 0);
    std::sort(v.begin(), v.end(),
              [vec] (size_t i, size_t j) {
        return vec[i] > vec[j];
    });

    return v;
}

}}  // namespace ant::std_ext
