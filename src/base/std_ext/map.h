#pragma once

#include <algorithm>
#include <type_traits>

#include "container.h"

namespace ant {
namespace std_ext {


// templates to match map-like containers
template <typename T>
struct is_pair : std::false_type {};

template <typename T, typename U>
struct is_pair<std::pair<T, U>> : std::true_type {};

template <typename...>
struct is_mapping : std::false_type {};

template <typename Container>
struct is_mapping<Container, std::enable_if<
    is_pair<typename std::iterator_traits<typename Container::iterator>::value_type>::value
>> : std::true_type {};



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
template <typename Map, typename = std::enable_if<is_mapping<Map>::value>>
typename Map::iterator min_map_element(Map& m)
{
    return std::min_element(m.begin(), m.end(), [] (typename Map::value_type& l,
                                                    typename Map::value_type& r) -> bool {
        return l.second < r.second;
    });
}

/**
 * @brief Takes map-like container and returns an iterator to the pair with the maximum associated key value
 * @param map map-like container with comparable objects
 * @return map::iterator to element with max value
 */
template <typename Map, typename = std::enable_if<is_mapping<Map>::value>>
typename Map::iterator max_map_element(Map& m)
{
    return std::max_element(m.begin(), m.end(), [] (typename Map::value_type& l,
                                                    typename Map::value_type& r) -> bool {
        return l.second < r.second;
    });
}

}}  // namespace ant::std_ext
