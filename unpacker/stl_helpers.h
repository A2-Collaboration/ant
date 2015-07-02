#ifndef STL_HELPERS_H
#define STL_HELPERS_H

#include <string>

namespace std_ext {

inline bool string_ends_with(std::string const& value, std::string const& ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

}

#endif // STL_HELPERS_H