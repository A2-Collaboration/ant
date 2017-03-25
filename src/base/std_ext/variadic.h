#pragma once

#include <cstddef> // for std::size_t
#include <tuple>

namespace ant {
namespace std_ext {

// simple index_sequence, as it's not available in C++11
template <std::size_t... Is>
struct indices {};

template <std::size_t N, std::size_t... Is>
struct build_indices : build_indices<N-1, N-1, Is...> {};

template <std::size_t... Is>
struct build_indices<0, Is...> : indices<Is...> {};

// test if T is a specialization of template
template <typename T, template <typename...> class Template>
struct is_specialization_of : std::false_type {};
template <template <typename...> class Template, typename... Args>
struct is_specialization_of<Template<Args...>, Template> : std::true_type {};

// helper for tuple
template<typename T, typename... Ts, std::size_t... Idx>
constexpr
std::tuple<Ts...> strip_first_from_tuple_impl(const std::tuple<T, Ts...>& t, indices<Idx...>) noexcept {
    return std::make_tuple(std::get<Idx+1>(t)...);
}

template<typename T, typename... Ts>
constexpr
std::tuple<Ts...> strip_first_from_tuple(const std::tuple<T, Ts...>& t) noexcept {
    return strip_first_from_tuple_impl(t, build_indices<sizeof...(Ts)>());
}

// extent std::tuple_element to return void if not present
template<std::size_t I, typename Tuple>
struct tuple_element_void {
    using type = typename std::tuple_element<I, Tuple>::type;
};

template<std::size_t I>
struct tuple_element_void<I, std::tuple<>> {
    using type = void;
};


template<typename Arg = void, typename... Args>
struct first_element_is_tuple {
    static constexpr bool value = std_ext::is_specialization_of<Arg, std::tuple>::value;
};

}}
