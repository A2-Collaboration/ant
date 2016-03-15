#pragma once

#include <cmath>

// ROOT compat
#include "TVector3.h"

namespace ant {

struct vec3 {
    double x = {};
    double y = {};
    double z = {};

    vec3() noexcept = default;
    vec3(const vec3&) noexcept = default;
    vec3(vec3&&) noexcept  = default;
    vec3(const double X, const double Y, const double Z) noexcept : x(X), y(Y), z(Z) {}

    vec3& operator= (const vec3&) noexcept = default;
    vec3& operator= (vec3&&) noexcept = default;


    // =====  TVector3 interface =====

    vec3(const TVector3& v) noexcept:
        x(v.X()),y(v.Y()),z(v.Z())
    {}

    vec3& operator= (const TVector3& v) noexcept {
        x = v.X();
        y = v.Y();
        z = v.Z();
        return *this;
    }

    operator TVector3() const {
        return TVector3(x,y,z);
    }

    // ===============================

    vec3& operator+=(const vec3& other) noexcept {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    vec3& operator-=(const vec3& other) noexcept {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    vec3 operator+(const vec3& other) const noexcept {
        return vec3(*this) += other;
    }

    vec3 operator-(const vec3& other) const noexcept {
        return vec3(*this) -= other;
    }

    vec3& operator*=(const double s) noexcept {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    vec3& operator/=(const double s) noexcept {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    vec3 operator*(const double s) const noexcept {
        return vec3(*this) *= s;
    }

    vec3 operator/(const double s) const noexcept {
        return vec3(*this) /= s;
    }

    double R2() const noexcept {
        return x*x+y*y+z*z;
    }

    double R() const noexcept {
        return sqrt(this->R2());
    }

    /**
     * @brief Angle
     * @param other
     * @return
     *
     * Taken from ROOT TVector3
     */
    double Angle(const vec3& other) const {

        const auto ptot2 = this->R2() * other.R2();

        if(ptot2 <= 0) {
           return 0.0;
        } else {
           auto arg = this->Dot(other) / sqrt(ptot2);
           if(arg >  1.0) arg =  1.0;
           if(arg < -1.0) arg = -1.0;
           return acos(arg);
        }
    }

    double Dot(const vec3& other) const noexcept {
        return x * other.x + y * other.y + z * other.z;
    }

    vec3 Cross(const vec3& p) const noexcept {
       return vec3(y*p.z-p.y*z, z*p.x-p.z*x, x*p.y-p.x*y);
    }

    double Phi() const {
        return x == 0.0 && y == 0.0 ? 0.0 : atan2(y,x);
    }

    /**
     * @brief Theta
     * @return
     * Taken from TVector3
     */
    double Theta() const {
        return x == 0.0 && y == 0.0 && z == 0.0 ? 0.0 : atan2(sqrt(x*x+y*y),z);
    }

    vec3 Unit() const
    {
       // return unit vector parallel to this.
       const auto tot2 = this->R2();
       const auto tot = (tot2 > 0) ?  1.0/sqrt(tot2) : 1.0;
       return vec3(*this)*=tot;
    }

    static vec3 RThetaPhi(const double r, const double theta, const double phi) {
        const auto stheta = sin(theta);
        const auto ctheta = cos(theta);
        const auto sphi   = sin(phi);
        const auto cphi   = cos(phi);
        return vec3(r*stheta*cphi, r*stheta*sphi, r*ctheta);
    }

    bool operator==(const vec3& other) const noexcept {
        return x==other.x && y == other.y && z == other.z;
    }

    bool operator!=(const vec3& other) const noexcept {
        return !(*this == other);
    }

};

}
