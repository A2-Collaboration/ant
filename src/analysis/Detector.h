#pragma once

#include "base/printable.h"

namespace ant {

class detector_t : public printable_traits {
private:
    unsigned int v;

public:
    detector_t() : v(0) {}
    detector_t(unsigned int _v) : v(_v) {}

    const static detector_t None;
    const static detector_t CB;
    const static detector_t PID;
    const static detector_t MWPC0;
    const static detector_t MWPC1;
    const static detector_t TAPS;
    const static detector_t TAPSVeto;
    const static detector_t Cherenkov;

    const static detector_t MWPC;
    const static detector_t anyCB;
    const static detector_t anyTAPS;
    const static detector_t anyVeto;

    bool operator==(const detector_t& o) const {
        return v==o.v;
    }

    bool operator!=(const detector_t& o) const {
        return v!=o.v;
    }

    operator bool() const {
        return v;
    }

    detector_t operator&(const detector_t& o) const {
        return v & o.v;
    }

    detector_t operator^(const detector_t& o) const {
        return v ^ o.v;
    }

    detector_t operator^=(const detector_t&o) {
        v ^= o.v;
        return *this;
    }

    bool operator<(const detector_t& o) const {
        return v < o.v;
    }

    bool operator>(const detector_t& o) const {
        return v > o.v;
    }

    detector_t& operator|=(const detector_t& o) {
        v |= o.v;
        return *this;
    }

    detector_t operator|(const detector_t& o) const {
        return  detector_t(v | o.v);
    }

    std::ostream& Print(std::ostream &stream) const override;

    operator std::string() const;
};

}
