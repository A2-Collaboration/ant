#pragma once

#include "interval.h"
#include "base/std_ext/printable.h"

#include <vector>

namespace ant {

template <typename T>
class PiecewiseInterval : public std::vector<interval<T>>
{

public:

    using interval_t = ant::interval<T>;

    // ctor
    using std::vector<interval_t>::vector;

    /**
     * @brief Test if any interval contains x
     * @param x Value to check
     * @return contained yes/no
     */
    bool Contains(const T& x) const {
        for(auto& i : *this) {
            if(i.Contains(x))
                return true;
        }
        return false;
    }

    /**
     * @brief Simplify the list by joining intervals that overlap
     * May reduce the number of intervals in the list and speed up following Contains() checks a bit
     */
    void Compact() {
        auto i =this->begin();
        while(i!= this->end()) {
            auto j = i;
            ++j;
            while(j!=this->end()) {

                auto res = i->tryJoinWith(*j);

                if( res )
                    j = this->erase(j);
                else
                    ++j;
            }
            ++i;
        }
    }

    interval_t EnclosingInterval() const {
        if(this->empty())
            return interval_t({},{});

        auto min_start =  [] (const interval_t& a, const interval_t& b) {
            return a.Start() < b.Start();
        };
        auto min = std::min_element(this->begin(), this->end(), min_start);
        auto max_stop =  [] (const interval_t& a, const interval_t& b) {
            return a.Stop() < b.Stop();
        };
        auto max = std::max_element(this->begin(), this->end(), max_stop);

        return {min->Start(), max->Stop()};
    }

    /**
     * @brief Check equality.
     * Two PiecewiseIntervals are equal if they contain exactly the same intervals
     * @param rhs
     * @return
     */
    bool operator== (const PiecewiseInterval<T>& rhs) const {
        auto i = this->begin();
        auto j = rhs.begin();

        while(i!=this->end() && j != rhs.end()) {
            if( *i != *j )
                return false;
            ++i;
            ++j;
        }
        return i==this->end() && j==rhs.end();
    }

    friend std::ostream& operator<<(std::ostream& stream, const PiecewiseInterval& o) {
        stream << "[";
        for(const auto& i : o) {
            stream << i;
        }
        stream << "]";
        return stream;
    }

    friend std::istream& operator>>(std::istream& out, PiecewiseInterval<T>& t)
    {
        out >> std::ws;
        if(out.peek() == '[')
            out.ignore();
        t.clear();

        interval<T> in({}, {});
        while(out.peek() != EOF && out >> in) {
            t.push_back(in);
            if(out.peek() == ';')
                out.ignore();
        }
        if(!t.empty())
            out.clear();

        return out;
    }

    T Area() const {
        T sum = 0.0;
        for(const auto& i: *this) {
            sum += i.Length();
        }
        return sum;
    }
};

// prevent matching as printable container
namespace std_ext {
template<typename T>
struct is_stl_container_like<ant::PiecewiseInterval<T>>
{
    static constexpr bool value = false;
};
}

} // namespace ant
