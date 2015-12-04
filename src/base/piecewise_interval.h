#pragma once

#include "interval.h"
#include "printable.h"

#include <vector>

namespace ant {

template <typename T>
class PiecewiseInterval : public std::vector<interval<T>>, public printable_traits
{

public:

    using std::vector<ant::interval<T>>::vector;

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

    std::ostream& Print(std::ostream& stream) const {
        stream << "[";
        for(const auto& i : *this) {
            stream << i;
        }
        stream << "]";
        return stream;
    }

    friend std::istream& operator>>(std::istream& out, PiecewiseInterval<T>& t)
    {
        // try

        out >> std::ws;
        if(out.peek() == '[')
            out.ignore();
        t.resize(0);

        interval<T> in;
        while(out.peek() != EOF && out >> in) {
            t.push_back(in);
            if(out.peek() == ';')
                out.ignore();
        }
        if(!t.empty())
            out.clear();

        return out;
    }
};

}
