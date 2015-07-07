#include "printable.h"

std::ostream& operator<< (std::ostream& stream, const ant::printable_traits* printable) {
    return printable->Print(stream);
}

std::ostream& operator<< (std::ostream& stream, const ant::printable_traits& printable) {
    return printable.Print(stream);
}
