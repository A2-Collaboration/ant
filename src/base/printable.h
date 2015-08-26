#pragma once

#include <vector>
#include <ostream>

namespace ant {

class printable_traits {
public:
  virtual std::ostream& Print( std::ostream& stream ) const =0;
  virtual ~printable_traits() = default;
};

}

inline std::ostream& operator<< (std::ostream& stream, const ant::printable_traits& printable) {
    return printable.Print(stream);
}

template<class T>
inline std::ostream& operator<< (std::ostream& stream, const std::vector<T>& v)
{
    for (auto& entry: v)
    {
        stream << entry << " , ";
    }
    return stream;
}
