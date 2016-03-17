#include "LorentzVec.h"

#include "TLorentzVector.h"

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
