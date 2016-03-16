#pragma once

#include "vec2.h"

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

    /**
     * @brief Get the length squared
     * @return
     * @see R()
     */
    double R2() const noexcept {
        return x*x+y*y+z*z;
    }

    /**
     * @brief Get the length or radius
     * @return
     * @see R2()
     */
    double R() const noexcept {
        return sqrt(this->R2());
    }

    /**
     * @brief Angle to another vec3
     * @param other
     * @return (radians)
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

    /**
     * @brief Dot product
     * @param other
     * @return
     */
    double Dot(const vec3& other) const noexcept {
        return x * other.x + y * other.y + z * other.z;
    }

    /**
     * @brief Cross product
     * @param p other vector
     * @return
     */
    vec3 Cross(const vec3& p) const noexcept {
       return vec3(y*p.z-p.y*z, z*p.x-p.z*x, x*p.y-p.x*y);
    }

    /**
     * @brief Calculate the phi angle
     * @return (radians)
     */
    double Phi() const {
        return x == 0.0 && y == 0.0 ? 0.0 : atan2(y,x);
    }

    /**
     * @brief Calculate the theta angle
     * @return (radians)
     * Taken from TVector3
     */
    double Theta() const {
        return x == 0.0 && y == 0.0 && z == 0.0 ? 0.0 : atan2(sqrt(x*x+y*y),z);
    }

    /**
     * @brief Get a vector of unit length parallel to this one
     * @return
     */
    vec3 Unit() const noexcept
    {
       // return unit vector parallel to this.
       const auto tot2 = this->R2();
       const auto tot = (tot2 > 0) ?  1.0/sqrt(tot2) : 1.0;
       return vec3(*this)*=tot;
    }

    /**
     * @brief Create a new vec3 from Radius and theta and phi angles
     * @param r Radius/Length of vector
     * @param theta (radians) Angle to the z-axis
     * @param phi (radians)
     * @return
     */
    static vec3 RThetaPhi(const double r, const double theta, const double phi) noexcept {
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

    /**
     * @brief Get X and Y as vec2
     * @return (X,Y)
     */
    vec2 XY() const noexcept {
        return {x,y};
    }

    /**
     * @brief Get Y and Z as a vec2
     * @return (Y,Z)
     */
    vec2 YZ() const noexcept {
        return {y,z};
    }

    /**
     * @brief Get X and Z as a vec2
     * @return (X,Z)
     */
    vec2 XZ() const noexcept {
        return {x,z};
    }

};

}
