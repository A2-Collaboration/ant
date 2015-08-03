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

std::ostream& operator<< (std::ostream& stream, const ant::printable_traits& printable);

template<class T>
std::ostream& operator<< (std::ostream& stream, const std::vector<T>& v)
{
    for (auto& entry: v)
    {
        stream << entry << " , ";
    }
    return stream;
}

//std::ostream& operator<< (std::ostream& stream, const ant::printable_traits* printable);
