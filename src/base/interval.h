#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <algorithm>
#include <stdexcept>
#include "printable.h"
#include "types.h"

namespace ant {

/**
 * @brief Interval class template
 */
template <class T> class interval: public printable_traits {
protected:
    T _start;
    T _stop;

public:
    interval( const T start=0.0, const T stop=0.0 ): _start(start), _stop(stop) {}

    /**
     * @brief Factory Function: Create interval from center and widht
     * @param center the center
     * @param width the width
     * @return an interval that spans symmetrically around center with given width
     */
    static interval CenterWidth( const T center, const T width) {
        return std::move( interval(center - width/2.0, center + width/2.0));
    }

    virtual ~interval() {}

    /**
     * @brief Get the lower boundary
     * @return start position
     */
    const T& Start() const { return _start; }

    /**
     * @brief Get the lower boundary
     * @return start position
     */
    T& Start() { return _start; }

    /**
     * @brief Get the upper boundary
     * @return
     */
    const T& Stop()  const { return _stop; }

    /**
     * @brief Get the upper boundary
     * @return
     */
    T& Stop() { return _stop; }

    /**
     * @brief Get the length of the interval
     * @return length
     */
    T Length() const { return Stop() - Start(); }

    /**
     * @brief Get the center
     * @return center
     */
    T Center() const { return (Start() + Stop()) / 2.0; }

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

    bool Contains( const T& x ) const { return _start <= x && _stop >= x; }

    bool Disjoint( const ant::interval<T>& i ) const {
        return (Stop() < i.Start()) || (i.Stop() < Start());
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

    bool IsSane() const { return _start <= _stop; }

    // printable_traits interface
    std::ostream &Print(std::ostream &stream) const {
        stream << "[" << _start << ":" << _stop << "]";
        return stream;
    }

};

/**
 * @brief Interval of doubles
 */
typedef ant::interval<double> IntervalD;

/**
 * @brief Interval of ints
 */
typedef ant::interval<int> IntervalI;
}

#endif /* INTERVAL_H_ */
