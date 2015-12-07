#include "Event.h"
#include <iostream>

using namespace std;
using namespace ant::analysis::data;

const ParticleList Event::Data::PTypeLists::empty;


std::ostream &Event::Print(std::ostream &stream) const
{

    stream << "----- Event ---------------------------------\n";

    stream << "== Reconstructed ==\n";
    stream << Reconstructed << '\n';
    stream << "== MC True ==\n";
    stream << MCTrue << '\n';

    stream << "---------------------------------------------\n";

    return stream;

}


std::ostream &Event::Data::Print(std::ostream &stream) const
{
    stream << "-- Particles --\n";
    stream << Particles << "\n";
    stream << "-- Intermediates --\n";
    stream << Intermediates << "\n";

    stream << "-- Candidates --\n";
    for( auto& candidate : Candidates ) {
        stream << *candidate << '\n';
    }

    stream << "-- Tagger --\n";
    for( auto& taggerhit : TaggerHits ) {
        stream << taggerhit << '\n';
    }

    stream << "-- Trigger --\n";
    cout << Trigger << "\n";

    return stream;

}


std::ostream &Event::Data::PTypeLists::Print(std::ostream &stream) const
{
    for( auto& particle : particles ) {
        stream << *particle << '\n';
    }

    return stream;
}
