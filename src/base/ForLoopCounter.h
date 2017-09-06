#pragma once
#include <iterator>
#include <iostream>
#include <algorithm>

namespace ant {

template <typename T>
struct ForLoopCounter {
    const T start;
    const T stop;

    constexpr ForLoopCounter(T Stop) noexcept : start(0), stop(Stop) {}
    constexpr ForLoopCounter(T Start, T Stop) noexcept : start(Start), stop(Stop) {}

    constexpr T counts() const noexcept {
        return stop-start;
    }

    struct iterator : std::iterator_traits<T> {
        const ForLoopCounter& counter;
        T i;
        T last;
        iterator(const ForLoopCounter& c, T start) noexcept : counter(c), i(start), last(i) {}

        iterator& operator++() {
            ++i;
            if(i-last >= 0.1*counter.counts()) {
                std::cout << "Progress: " << 100.0*(i-counter.start)/counter.counts() << "%" << std::endl;
                last = i;
            }
            return *this;
        }

        bool constexpr operator==(const iterator& other) const noexcept {
            return i == other.i;
        }

        bool constexpr operator!=(const iterator& other) const noexcept {
            return i != other.i;
        }

        const constexpr T& operator*() const noexcept {
            return i;
        }

    };

    constexpr iterator begin() const noexcept {
        return {*this, start};
    }

    constexpr iterator end() const noexcept {
        return {*this, stop};
    }

    constexpr iterator cbegin() const noexcept {
        return {*this, start};
    }

    constexpr iterator cend() const noexcept {
        return {*this, stop};
    }
};

//template<typename T, typename F>
//void ForEachCounter(const T& t, const F& f)
//{
//    const ForLoopCounter<T> counter(t);

//    std::for_each(counter.begin(),counter.end(),
//                  f());
//}

}
