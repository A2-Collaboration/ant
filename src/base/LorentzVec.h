#pragma once

#include "base/vec3.h"

// ROOT compat
#include "TLorentzVector.h"


namespace ant {

/**
 * @brief A Lorentz vector (x,y,z,E), diag(-1,-1,-1, 1)
 */
struct LorentzVec {
    vec3   x = {};
    double E = {};

    LorentzVec() noexcept = default;
    LorentzVec(const LorentzVec&) noexcept = default;
    LorentzVec(LorentzVec&&) noexcept = default;
    LorentzVec(const vec3& v, const double E_) noexcept:
        x(v), E(E_) {}
    LorentzVec(const double X, const double Y, const double Z, const double _E) noexcept :
        x(X,Y,Z), E(_E) {}

    LorentzVec& operator=(const LorentzVec&) noexcept = default;
    LorentzVec& operator=(LorentzVec&&) noexcept = default;


    // ====== TLorentzVector interface ======

    operator TLorentzVector() const {
        return TLorentzVector(x, E);
    }

    LorentzVec(const TLorentzVector& other) noexcept:
        x(other.Vect()), E(other.E()) {}

    LorentzVec& operator=(const TLorentzVector other) noexcept {
        x = other.Vect();
        E = other.E();
        return *this;
    }

    // ======================================

    static LorentzVec EPThetaPhi(const double E, const double p, const double theta, const double phi) noexcept {
        return LorentzVec(vec3::RThetaPhi(p,theta,phi), E);
    }

    LorentzVec& operator+=(const LorentzVec& other) noexcept {
        x += other.x;
        E += other.E;
        return *this;
    }

    LorentzVec& operator-=(const LorentzVec& other) noexcept {
        x -= other.x;
        E -= other.E;
        return *this;
    }

    LorentzVec operator+(const LorentzVec& other) const noexcept {
        return LorentzVec(*this) += other;
    }

    LorentzVec operator-(const LorentzVec& other) const noexcept {
        return LorentzVec(*this) -= other;
    }

    double M2() const noexcept {
        return E*E - x.R2();
    }

    double M() const {
        const auto mm = this->M2();
        return mm < 0.0 ? -sqrt(-mm) : sqrt(mm);
    }

    double Theta() const {
        return x.Theta();
    }

    double Phi() const {
        return x.Phi();
    }

    double P() const {
        return x.R();
    }

    bool operator==(const LorentzVec& other) const noexcept {
        return x == other.x && E == other.E;
    }

    bool operator!=(const LorentzVec& other) const noexcept {
        return !(*this == other);
    }

    LorentzVec& operator*=(const double a) noexcept {
        x *= a;
        E *= a;
        return *this;
    }

    LorentzVec& operator/=(const double a) noexcept {
        x /= a;
        E /= a;
        return *this;
    }

    LorentzVec operator* (const double a) const noexcept {
        return LorentzVec(*this)*=a;
    }

    LorentzVec operator/ (const double a) const noexcept {
        return LorentzVec(*this)/=a;
    }

    double Beta() const noexcept {
        return x.R() / E;
    }

    /**
     * @brief Gamma = 1 / sqrt(1-beta^2)
     * @return
     */
    double Gamma() const noexcept {
        return 1.0 / sqrt(1 - (x.R2()/(E*E)));
    }

    double Dot(const LorentzVec& other) const noexcept {
        return E*other.E - x.Dot(other.x);
    }

    vec3 BoostVector() const noexcept {
        return x / E;
    }

};

}
