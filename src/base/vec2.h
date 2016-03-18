#pragma once

#include <cmath>

class TVector2;

namespace ant {

struct vec2 {
    double x = {};
    double y = {};

    vec2() noexcept = default;
    vec2(const vec2&) noexcept = default;
    vec2(vec2&&) noexcept  = default;
    vec2(const double X, const double Y) noexcept : x(X), y(Y) {}

    vec2& operator= (const vec2&) noexcept = default;
    vec2& operator= (vec2&&) noexcept = default;


    // =====  TVector2 interface =====

    vec2(const TVector2& v) noexcept;
    vec2& operator=(const TVector2& v) noexcept;
    operator TVector2() const noexcept;

    // ===============================

    vec2& operator+=(const vec2& other) noexcept {
        x += other.x;
        y += other.y;
        return *this;
    }

    vec2& operator-=(const vec2& other) noexcept {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    vec2 operator+(const vec2& other) const noexcept {
        return vec2(*this) += other;
    }

    vec2 operator-(const vec2& other) const noexcept {
        return vec2(*this) -= other;
    }

    vec2& operator*=(const double s) noexcept {
        x *= s;
        y *= s;
        return *this;
    }

    vec2& operator/=(const double s) noexcept {
        x /= s;
        y /= s;
        return *this;
    }

    vec2 operator*(const double s) const noexcept {
        return vec2(*this) *= s;
    }

    vec2 operator/(const double s) const noexcept {
        return vec2(*this) /= s;
    }

    /**
     * @brief Get the length or radius squared
     * @return
     * @see R()
     */
    double R2() const noexcept {
        return x*x+y*y;
    }

    /**
     * @brief  Get the length or Radius of the vector
     * @return length
     * @see R2()
     */
    double R() const noexcept {
        return sqrt(this->R2());
    }

    /**
     * @brief Calcualte the angle between two vec2.
     * @param other
     * @return angle in radians
     *
     * Taken from ROOT TVector3
     */
    double Angle(const vec2& other) const {

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
    double Dot(const vec2& other) const noexcept {
        return x * other.x + y * other.y;
    }

    /**
     * @brief Phi angle
     * @return radians
     */
    double Phi() const {
        return x == 0.0 && y == 0.0 ? 0.0 : atan2(y,x);
    }

    /**
     * @brief vector parallel to this one with R() == 1.0
     * @return unit length vector
     */
    vec2 Unit() const noexcept
    {
       // return unit vector parallel to this.
       const auto tot2 = this->R2();
       const auto tot = (tot2 > 0) ?  1.0/sqrt(tot2) : 1.0;
       return vec2(*this)*=tot;
    }

    /**
     * @brief Create vec2 from Radius and Angle
     * @param r Radius
     * @param phi Angle in radians
     * @return new vec2
     */
    static vec2 RPhi(const double r, const double phi) noexcept {
        const auto sphi   = sin(phi);
        const auto cphi   = cos(phi);
        return vec2(r*cphi, r*sphi);
    }

    bool operator==(const vec2& other) const noexcept {
        return x==other.x && y == other.y;
    }

    bool operator!=(const vec2& other) const noexcept {
        return !(*this == other);
    }

    /**
     * @brief Phi_mpi_pi returns phi angle in interval [-pi,pi)
     * @param phi
     * @return
     * @note copied from TVector2::Phi_mpi_pi
     */
    static double Phi_mpi_pi(double phi) {
       while (phi >= M_PI) phi -= 2*M_PI;
       while (phi < -M_PI) phi += 2*M_PI;
       return phi;
    }

};

}
