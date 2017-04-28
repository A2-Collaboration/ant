#pragma once

#include "types.h"

#include <algorithm>
#include <stdexcept>
#include <limits>
#include <istream>
#include <sstream>
#include <cmath>

namespace ant {

/**
 * @brief Interval class template
 */
template<class T>
class interval {
protected:
    T _start;
    T _stop;

public:
    constexpr interval(const T& start, const T& stop) noexcept :
        _start(start), _stop(stop)
    {}


    /**
     * @brief Factory Function: Create interval from center and width
     * @param center the center
     * @param width the width
     * @return an interval that spans symmetrically around center with given width
     */
    constexpr static interval CenterWidth( const T center, const T width) noexcept {
        return interval(center - width/2.0, center + width/2.0);
    }

    /**
     * @brief Get the lower boundary
     * @return start position
     */
    constexpr const T& Start() const noexcept { return _start; }

    /**
     * @brief Get the lower boundary
     * @return start position
     */
    T& Start() noexcept { return _start; }

    /**
     * @brief Get the upper boundary
     * @return
     */
    constexpr const T& Stop() const noexcept { return _stop; }

    /**
     * @brief Get the upper boundary
     * @return
     */
    T& Stop() noexcept { return _stop; }

    /**
     * @brief Get the length of the interval
     * @return length
     */
    constexpr T Length() const noexcept { return Stop() - Start(); }

    /**
     * @brief Get the center
     * @return center
     */
    constexpr T Center() const noexcept { return (Start() + Stop()) / 2.0; }

    /**
     * @brief Add a value to both boundaries. This shifts the interval to higher values.
     * @param a The value to add
     * @return itself (Interval)
     */
    interval<T>& operator+= ( const T& a ) { _start += a; _stop += a; return *this; }

    /**
     * @brief Subtract a value to both boundaries. This shifts the interval to lower values.
     * @param a The value to subtract
     * @return itself (Interval)
     */
    interval<T>& operator-= ( const T& a ) { _start -= a; _stop -= a; return *this; }

    /**
     * @brief Multiply both boundaries with a factor. This scales the interval (but also moves the center if center != 0)
     * @param f The factor to multiply with
     * @return itself (Interval)
     */
    interval<T>& operator*= ( const T& f ) { _start *= f; _stop *= f; return *this; }

    /**
     * @brief Deivide both boundaries by a factor. This scales the interval (but also moves the center if center != 0)
     * @param f The factor to divide by
     * @return itself (Interval)
     */
    interval<T>& operator/= ( const T& f ) { _start /= f; _stop /= f; return *this; }

    /**
     * @brief Add a value to both boundaries and return the resulting interval (a new copy)
     *  The result is shifted to higher values (if a positive, of course).
     * @param a The value to add
     * @return this + a
     */
    interval<T> operator+ ( const T& a ) const { interval<T> tmp(*this); return tmp += a; }

    /**
     * @brief Subtract a value to both boundaries and return the resulting interval (a new copy)
     *  The result is shifted to lower values (if a positive, of course).
     * @param a The value to subtract
     * @return this - a
     */
    interval<T> operator- ( const T& a ) const { interval<T> tmp(*this); return tmp -= a; }

    /**
     * @brief Multiply both boundaries with a factor and return the resulting interval (a new copy).
     *  The result is scaled and the center is shifted if != 0.
     * @param f The factor to multiply with
     * @return this * f
     */
    interval<T> operator* ( const T& f ) const { interval<T> tmp(*this); return tmp *= f; }

    /**
     * @brief Deivide both boundaries by a factor and return the resulting interval (a new copy).
     *  The result is scaled and the center is shifted if != 0.
     * @param f The factor to divide by
     * @return itself (Interval)
     */
    interval<T> operator/ ( const T& f ) const { interval<T> tmp(*this); return tmp /= f; }


    /**
     * @brief operator []: Access boundaries array-style
     * @param n =0 -> lower, =1 -> upper
     * @return boundary n
     * @throw std::out_of_range if index != {0,1}
     */
    T operator[] ( const index_t n ) const {
        switch (n) {
        case 0:
            return _start;
        case 1:
            return _stop;
        default:
            throw std::out_of_range ("Interval Index invalid");
        }
    }

    constexpr bool operator == (const interval<T>& rhs) const noexcept {
        return Start() == rhs.Start() && Stop() == rhs.Stop();
    }

    constexpr bool operator != (const interval<T>& rhs) const noexcept {
        return !(*this == rhs);
    }

    void SetWidth(const T& width) {
        const T center = Center();
        Start() = center - width/2.0;
        Stop() = center + width/2.0;
    }

    void SetCenter(const T& center) {
        const T width = Length();
        Start() = center - width/2.0;
        Stop() = center + width/2.0;
    }

    constexpr bool Contains( const T& x ) const noexcept { return _start <= x && x <= _stop; }

    constexpr bool Disjoint( const ant::interval<T>& i ) const noexcept {
        return (Stop() < i.Start()) || (i.Stop() < Start());
    }

    /**
     * @brief Extend the interval if the given covers a larger range
     * @param other another interval
     */
    void Extend(const interval<T>& other) {
        Start() = std::min(other.Start(), Start());
        Stop()  = std::max(other.Stop(),  Stop());
    }

    /**
     * @brief Extend the interval if the given covers a larger range
     * @param other another interval
     */
    void Extend(const T& other) {
        Start() = std::min(other, Start());
        Stop()  = std::max(other, Stop());
    }

    /**
     * @brief Join another interval into this one if they overlap
     * @param other another interval
     * @return true if joined, false if intervals are disjoint
     */
    bool tryJoinWith(const ant::interval<T>& other) {
        if( ! Disjoint(other) ) {
            Extend(other);
            return true;
        }
        return false;
    }


    /**
     * @brief Make sure, start < stop
     *
     * This protected method is called whenever a write access to _start or _stop occurs.
     */
    void MakeSane() {
        if( _stop < _start )
            std::swap( _start, _stop );
    }

    constexpr bool IsSane() const noexcept { return _start <= _stop; }

    static constexpr interval getMaxRange() noexcept {
        return std::numeric_limits<T>::has_infinity ?
                    interval(-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity() ) :
                    interval(std::numeric_limits<T>::min(), std::numeric_limits<T>::max() );
    }

    static constexpr interval getMaxPositiveRange() noexcept {
        return std::numeric_limits<T>::has_infinity ?
                    interval( 0, std::numeric_limits<T>::infinity() ) :
                    interval( 0, std::numeric_limits<T>::max() );
    }

    static constexpr interval getMaxNegativeRange() noexcept {
        return std::numeric_limits<T>::has_infinity ?
                    interval( -std::numeric_limits<T>::infinity(), 0 ) :
                    interval( std::numeric_limits<T>::min(), 0 );
    }


    // default stringification
    friend std::ostream& operator<<(std::ostream& stream, const interval<T>& i) {
        return stream << "[" << i._start << ":" << i._stop << "]";
    }

    // the >> operator parses stringified versions
    friend std::istream& operator>>(std::istream& out, interval<T>& t)
    {
        // skip leading whitespace
        out >> std::ws;
        // skip leading [
        bool have_leading_bracket = false;
        if(out.peek() == '[') {
            out.ignore();
            have_leading_bracket = true;
        }
        // read Start
        if(out >> t.Start()) {
            // maybe read stop?
            if(out.peek() == ':' || out.peek() == '-') {
                out.ignore();
                out >> t.Stop();
            }
            else {
                t.Stop() = t.Start();
            }

            out >> std::ws;
            auto nextchar = out.peek();
            // check closing ]
            if(nextchar == ']') {
                if(have_leading_bracket)
                    out.ignore();
                else
                    out.setstate(std::ios::failbit);
            }
            else {
                // did not find closing ]
                if(have_leading_bracket)
                    out.setstate(std::ios::failbit);
                // prevent , as next char
                else if(nextchar == ',')
                    out.setstate(std::ios::failbit);
                // do not mark EOF as error
                else if(nextchar == EOF)
                    out.clear();
            }
        }
        return out;
    }

    /**
     * @brief Clip x to be inside the interval
     * @param x
     * @return
     */
    constexpr T Clip(const T& x) const {
        return x > _stop ? _stop : (x < _start ? _start : x);
    }

    /**
     * @brief Round boundaries to integer using std::round
     * @return reference to rounded instance
     * @note interval may become not sane by this operation
     */
    interval<T>& Round() noexcept {
        _start = std::round(_start);
        _stop = std::round(_stop);
        return *this;
    }

    /**
     * @brief AsRangeString makes "mathematical" string expression of interval
     * @param label the label this interval applies to
     * @return the string
     * @note it works well with _start, _stop being infinite/nan
     */
    std::string AsRangeString(const std::string& label = "x") const noexcept {
        std::stringstream ss;
        if(std::isfinite(_start))
            ss << _start << "<=";
        ss << label;
        if(std::isfinite(_stop))
            ss << "<=" << _stop;
        return ss.str();
    }
};

/**
 * @brief Interval of doubles
 */
using IntervalD = interval<double>;

/**
 * @brief Interval of ints
 */
using IntervalI = interval<int>;

template <typename T>
interval<T> intersect(const interval<T>& a, const interval<T>& b) {
    if(!a.Disjoint(b)) {
        return ant::interval<T> (
                    std::max(a.Start(), b.Start()),
                    std::min(a.Stop(), b.Stop())
                    );
    }
    return interval<T>({}, {});
}


} // namespace ant


