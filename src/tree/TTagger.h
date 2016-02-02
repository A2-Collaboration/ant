#pragma once

#include "TTaggerHit.h"

namespace ant {

struct TTagger : printable_traits
{
    std::vector<TTaggerHit> Hits;
    std::vector< TKeyValue<double> > Energies;

    TTagger() {}
    virtual ~TTagger() {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Hits, Energies);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        for(auto& hit : Hits)
            s << hit;
        return s;
    }

};

}
