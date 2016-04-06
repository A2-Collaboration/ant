#pragma once

#include <base/interval.h>
#include <cmath>

namespace ant {

template <typename T>
T step(const ant::interval<T> i, const unsigned n, const unsigned p) {
    return i.Start() + i.Length() / (n-1) * std::min(n-1,p);
}

struct stepper {
    const ant::interval<double>& i;
    const unsigned n;
    unsigned p=0;
    double value;

    stepper(const interval<double>& interv, const unsigned steps): i(interv), n(steps), value(step(i,n,0)) {}

    bool Next() {

        if(p<n) {
            ++p;
            value = step(i,n,p);
            return true;
        } else
            return false;
    }

    bool Done() const {
        return p >= n;
    }

};


template <typename T>
struct Range_t {

    Range_t(const ant::interval<T>& i, const unsigned n) noexcept:
        m_i(i), m_n(n) {}

    struct iterator {

        using value_type        = T;
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const T*;
        using reference         = const T;

        bool operator==(const iterator& other) const noexcept {
            return m_p == other.m_p;
        }

        bool operator!=(const iterator& other) const noexcept {
            return ! operator==(other);
        }

        iterator& operator++() noexcept {
            ++m_p; return *this;
        }

        iterator& operator--() noexcept {
            --m_p; return *this;
        }

        value_type value() const noexcept {
            return step(m_range.m_i, m_range.m_n, m_p);
        }

        unsigned position() const noexcept {
            return m_p;
        }

        reference operator*() const noexcept {
            return value();
        }

        iterator(const Range_t& range, const unsigned pos=0):
            m_range(range), m_p(pos) {}

        protected:
            const Range_t& m_range;
            unsigned m_p = 0;

    };

    iterator begin() const noexcept { return iterator(*this, 0); }
    iterator end()   const noexcept { return iterator(*this, m_n); }

protected:
    const ant::interval<T>& m_i;
    const unsigned m_n;

};

template <typename T>
Range_t<T> Range(const ant::interval<T>& i, const unsigned n) noexcept {
    return Range_t<T>(i,n);
}

}
