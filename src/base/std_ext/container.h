#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

namespace ant {
namespace std_ext {

template<typename T>
struct is_stl_container_like
{
    template<typename A>
    static constexpr bool test(
            A * pt,
            const A* cpt = nullptr,
            decltype(pt->begin()) * = nullptr,
            decltype(pt->end()) * = nullptr,
            decltype(cpt->begin()) * = nullptr,
            decltype(cpt->end()) * = nullptr)
    {

        typedef typename A::iterator iterator;
        typedef typename A::const_iterator const_iterator;
        return  std::is_same<decltype(pt->begin()),iterator>::value &&
                std::is_same<decltype(pt->end()),iterator>::value &&
                std::is_same<decltype(cpt->begin()),const_iterator>::value &&
                std::is_same<decltype(cpt->end()),const_iterator>::value;
    }

    template<typename A>
    static constexpr bool test(...) {
        return false;
    }

    typedef typename std::remove_const<T>::type test_type;
    static constexpr bool value = test<test_type>(nullptr);

};

// typical specialization to prevent matching std::string as a container
template<>
struct is_stl_container_like<std::string>
{
    static constexpr bool value = false;
};

template<typename Cont>
typename std::enable_if<is_stl_container_like<Cont>::value, void>::type
insertRange(Cont& v, unsigned start, unsigned stop) {
    int length = stop-start+1;
    if(length<1)
        return;
    typename std::decay<Cont>::type v_(static_cast<size_t>(length));
    std::iota(v_.begin(), v_.end(), start);
    v.insert(v.end(), v_.cbegin(), v_.cend());
}

template<typename Cont, typename U>
typename std::enable_if<is_stl_container_like<Cont>::value, void>::type
concatenate(Cont& dest, const U& src) {
    dest.insert(dest.end(), src.begin(), src.end());
}

template<typename Cont, typename U>
typename std::enable_if<is_stl_container_like<Cont>::value, bool>::type
contains(const Cont& v, const U& val) {
    return std::find(v.cbegin(), v.cend(), val) != v.cend();
}

}} // namespace ant::std_ext
