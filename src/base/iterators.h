#pragma once
#include <iterator>



template <typename T>
class circit: public std::iterator_traits<typename T::value_type> {
protected:
     const T start;
     const T end;
     T pos;
public:
    constexpr circit(T Start, T Stop) noexcept : start(Start), end(Stop), pos(start)  {}

    void next() {
        ++pos;
        if(pos==end)
            pos=start;
    }

    void previous() {
        if(pos==start) {
            pos=end;
        }
        --pos;
    }

    void Reset() {
        pos = start;
    }

    circit& operator++ () {
        next();
        return *this;
    }

    circit operator++(int) {circit tmp(*this); operator++(); return tmp; }

    constexpr typename T::value_type operator *() { return *pos; }

};

template <typename T>
circit<T> getCirculatIterator(T start, T end) {
    return circit<T>(start, end);
}
