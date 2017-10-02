#pragma once

#include <algorithm>
#include <type_traits>

#include "container.h"

namespace ant {
namespace std_ext {


// templates to match pair and map-like containers
template <typename T>
struct is_pair {
    static constexpr bool value = false;
};

template <typename T, typename U>
struct is_pair<std::pair<T, U>> {
    static constexpr bool value = true;
};

template <typename...>
struct is_mapping {
    static constexpr bool value = false;
};

template <typename Container>
struct is_mapping<Container, typename std::enable_if<
    is_pair<typename std::iterator_traits<typename Container::iterator>::second_type>::value
>> {
    static constexpr bool value = true;
};


// methods to access values of second_type from pairs, container of pairs, ...
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


template<class Pair>
typename std::enable_if<is_pair<Pair>::value, Pair>::type::second_type
get_second(Pair p) {
    return p.second;
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
        return get_second(l) < get_second(r);
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
        return get_second(l) < get_second(r);
    });
}

}}  // namespace ant::std_ext
