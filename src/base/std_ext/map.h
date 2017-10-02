#pragma once

#include <algorithm>
#include <type_traits>

#include "container.h"

namespace ant {
namespace std_ext {


// templates to match map-like containers
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
    is_pair<typename std::iterator_traits<typename Container::iterator>::value_type>::value
>> {
    static constexpr bool value = true;
};

template<class T, class = void>
struct has_second_type {
    static constexpr auto value = false;
};

template<class... Args> using void_t = typename voider<Args...>::type;

template<class T>
struct has_second_type<T, void_t<typename std::iterator_traits<typename T::iterator>::second_type>> {
    static constexpr auto value = true;
};


template<typename T>
struct pair_traits{
    using pair_type  = typename T::value_type;
    using key_type   = typename std::remove_const<typename pair_type::first_type>::type;
    using value_type = typename pair_type::second_type;
};

template<typename T>
using map_val_t = typename pair_traits<T>::value_type;



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


template <typename T>
using pair_t = typename pair_traits<T>::pair_type;

/**
 * @brief Takes map-like container and returns an iterator to the pair with the minimum associated key value
 * @param map map-like container with comparable objects
 * @return map::iterator to element with min value
 */
template <typename Map, typename = std::enable_if<has_second_type<Map>::value>>
typename Map::iterator min_map_element(Map& m)
{
    return std::min_element(m.begin(), m.end(), [] (pair_t<Map>& l,
                                                    pair_t<Map>& r) -> bool {
        return l.second < r.second;
    });
}

/**
 * @brief Takes map-like container and returns an iterator to the pair with the maximum associated key value
 * @param map map-like container with comparable objects
 * @return map::iterator to element with max value
 */
template <typename Map, typename = std::enable_if<has_second_type<Map>::value>>
typename Map::iterator max_map_element(Map& m)
{
    return std::max_element(m.begin(), m.end(), [] (pair_t<Map>& l,
                                                    pair_t<Map>& r) -> bool {
        return l.second < r.second;
    });
}

}}  // namespace ant::std_ext
