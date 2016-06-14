#include "LorentzVec.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include "TLorentzVector.h"
#pragma GCC diagnostic pop

using namespace ant;

LorentzVec::operator TLorentzVector() const noexcept {
    return TLorentzVector(p, E);
}

LorentzVec::LorentzVec(const TLorentzVector& other) noexcept :
    p(other.Vect()), E(other.E()) {}

LorentzVec& LorentzVec::operator=(const TLorentzVector other) noexcept {
    p = other.Vect();
    E = other.E();
    return *this;
}


