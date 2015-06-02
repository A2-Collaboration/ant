#include "Event.h"
#include <iostream>

using namespace  std;

const ant::ParticleList ant::Event::Data::PTypeLists::empty;


std::ostream &ant::Event::Print(std::ostream &stream) const
{

    stream << "------------------------------\n";

    stream << reconstructed << '\n';
    stream << mctrue << '\n';

    stream << "------------------------------\n";

    return stream;

}


std::ostream &ant::Event::Data::Print(std::ostream &stream) const
{
    stream << particles << "\n";
    stream << intermediates << "\n";

    for( auto& track : tracks ) {
        stream << track << '\n';
    }

    for( auto& taggerhit : taggerhits ) {
        stream << taggerhit << '\n';
    }

    cout << triggerinfo << "\n:";

    return stream;

}


std::ostream &ant::Event::Data::PTypeLists::Print(std::ostream &stream) const
{
    for( auto& particle : particles ) {
        stream << particle << '\n';
    }

    return stream;
}
