#ifndef PRINTABLE_H
#define PRINTABLE_H

#include <ostream>

namespace ant {

class printable_traits {
public:
    virtual std::ostream& Print( std::ostream& stream ) const =0;
};

}

std::ostream& operator<< (std::ostream& stream, const ant::printable_traits& printable);

std::ostream& operator<< (std::ostream& stream, const ant::printable_traits* printable);

#endif
