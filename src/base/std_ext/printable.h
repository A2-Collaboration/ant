#pragma once

#include "container.h"

#include <type_traits>
#include <ostream>

// put this into std namespace, as we target T exactly to
// those STL containers for ostream'ing (this is probably the only occasion
// where it makes sense to pollute std namespace)
// clang 3.9 otherwise complains when compiling test, btw
// (gcc 6.1 is happy though even if it's in namespace ant)
namespace std {

template<class T>
// use SFINAE to restrict this templated operator to STL containers such as vector,list,map,set
typename std::enable_if<ant::std_ext::is_stl_container_like<T>::value, std::ostream&>::type
operator<< (std::ostream& stream, const T& v)
{
    stream << "[";
    for(auto it = std::begin(v); it != std::end(v); ++it)
    {
        stream << *it;
        if(std::next(it) != std::end(v))
            stream << ";";
    }
    stream << "]";
    return stream;
}

// make std::pair printable for std::map support
template<class U, class V>
std::ostream&
operator<< (std::ostream& stream, const std::pair<U, V>& p) {
    return stream << p.first << "=" << p.second;
}

} // namespace std


