#include "Detector.h"
#include <ostream>
#include <sstream>

#include "expconfig/Detector_t.h"

using namespace std;
using namespace ant;


const detector_t detector_t::None(0);
const detector_t detector_t::CB(        Detector_t::ToBitfield(Detector_t::Type_t::CB));
const detector_t detector_t::PID(       Detector_t::ToBitfield(Detector_t::Type_t::PID));
const detector_t detector_t::MWPC0(     Detector_t::ToBitfield(Detector_t::Type_t::MWPC0));
const detector_t detector_t::MWPC1(     Detector_t::ToBitfield(Detector_t::Type_t::MWPC1));
const detector_t detector_t::TAPS(      Detector_t::ToBitfield(Detector_t::Type_t::TAPS));
const detector_t detector_t::TAPSVeto(  Detector_t::ToBitfield(Detector_t::Type_t::TAPSVeto));
const detector_t detector_t::Cherenkov( Detector_t::ToBitfield(Detector_t::Type_t::Cherenkov));

const detector_t detector_t::MWPC(    detector_t::MWPC0 | detector_t::MWPC1);
const detector_t detector_t::anyCB(   detector_t::CB    | detector_t::PID      | detector_t::MWPC);
const detector_t detector_t::anyTAPS( detector_t::TAPS  | detector_t::TAPSVeto);
const detector_t detector_t::anyVeto( detector_t::PID   | detector_t::TAPSVeto);

ostream& detector_t::Print(ostream& stream) const
{
    std::uint8_t i = 0;
    unsigned int c = 1;
    while(c<=v) {

        if(v & c) {

            stream << Detector_t::ToString(
                          static_cast<Detector_t::Type_t>(i)
                          )
                   << " ";
        }
        ++i;
        c <<= 1;
    }

    return stream;
}

ant::detector_t::operator std::string() const
{
    stringstream s;
    Print(s);
    return s.str();
}
