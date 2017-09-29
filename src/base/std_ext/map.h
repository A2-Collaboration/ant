#pragma once

#include <algorithm>
#include <type_traits>

#include "container.h"

namespace ant {
namespace std_ext {



template <typename Pair>
struct second_t {
    typename Pair::second_type operator()(const Pair& p) const {
        return p.second;
    }
};

template<typename Map>
second_t<typename Map::value_type>
second(const Map&) {
    return second_t<typename Map::value_type>();
}


/**
 * @brief Takes map-like container and returns an iterator to the pair with the minimum associated key value
 * @param map map-like container with comparable objects
 * @return map::iterator to element with min value
 */
template <typename T, typename = std::enable_if<std_ext::is_mapping<T>::value>>
typename T::iterator min_map_element(T& map)
{
    return std::min_element(map.begin(), map.end(), [] (typename T::value_type& l,
                                                        typename T::value_type& r) -> bool {
        return l.second < r.second;
    });
}

/**
 * @brief Takes map-like container and returns an iterator to the pair with the maximum associated key value
 * @param map map-like container with comparable objects
 * @return map::iterator to element with max value
 */
template <typename T, typename = std::enable_if<std_ext::is_mapping<T>::value>>
typename T::iterator max_map_element(T& map)
{
    return std::max_element(map.begin(), map.end(), [] (typename T::value_type& l,
                                                        typename T::value_type& r) -> bool {
        return l.second < r.second;
    });
}

}}  // namespace ant::std_ext
