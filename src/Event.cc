#include "Event.h"
#include <iostream>

using namespace  std;

const ant::ParticleList ant::Event::Data::PTypeLists::empty;


std::ostream &ant::Event::Print(std::ostream &stream) const
{

    stream << "----- Event ---------------------------------\n";

    stream << "== Reconstructed ==\n";
    stream << reconstructed << '\n';
    stream << "== MC True ==\n";
    stream << mctrue << '\n';

    stream << "---------------------------------------------\n";

    return stream;

}


std::ostream &ant::Event::Data::Print(std::ostream &stream) const
{
    stream << "-- Particles --\n";
    stream << particles << "\n";
    stream << "-- Intermediates --\n";
    stream << intermediates << "\n";

    stream << "-- Tracks --\n";
    for( auto& track : tracks ) {
        stream << *track << '\n';
    }

    stream << "-- Tagger --\n";
    for( auto& taggerhit : taggerhits ) {
        stream << *taggerhit << '\n';
    }

    stream << "-- Trigger --\n";
    cout << triggerinfo << "\n";

    return stream;

}


std::ostream &ant::Event::Data::PTypeLists::Print(std::ostream &stream) const
{
    for( auto& particle : particles ) {
        stream << *particle << '\n';
    }

    return stream;
}
