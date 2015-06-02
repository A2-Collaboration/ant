#include "Detector.h"
#include <ostream>
#include <sstream>

using namespace std;
using namespace ant;

const detector_t detector_t::None;
const detector_t detector_t::NaI(1);
const detector_t detector_t::PID(2);
const detector_t detector_t::MWPC(4);
const detector_t detector_t::BaF2(8);
const detector_t detector_t::PbWO4(16);
const detector_t detector_t::Veto(32);
const detector_t detector_t::anyCB(detector_t::NaI | detector_t::PID | detector_t::MWPC);
const detector_t detector_t::anyTAPS(detector_t::BaF2 | detector_t::PbWO4 | detector_t::Veto);
const detector_t detector_t::anyVeto(detector_t::PID | detector_t::Veto);

ostream& detector_t::Print(ostream& stream) const
{
    stream << "( ";
    if(*this & NaI) {
        stream << "NaI ";
    }
    if(*this & PID) {
        stream << "PID ";
    }
    if(*this & MWPC) {
        stream << "MWPC ";
    }
    if(*this & BaF2) {
        stream << "PbWO4 ";
    }
    if(*this & Veto) {
        stream << "Veto ";
    }
    stream << ")";
    return stream;
}

ant::detector_t::operator std::string() const
{
    stringstream s;
    Print(s);
    return s.str();
}
