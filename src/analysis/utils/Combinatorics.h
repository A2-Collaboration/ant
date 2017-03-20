#pragma once

#include "base/std_ext/container.h"

#include <vector>
#include <cmath>

namespace ant {
namespace analysis {
namespace utils {

/**
 * @brief KofNvector class: Generate all combinations when drawing k elements out of a vector.
 *
 * Elements are kept in the order they were in inside the data vector.
 *
 * No duplicate combinations.
 *
 * No duplicate elements (drawing without putting back).
 */
template <class T>
class NchooseK {

protected:
    const std::vector<T> _data;

    using index_type = std::size_t;
    using index_list = std::vector<index_type>;
    index_list indices;

    bool _done;


    bool nextlevel( index_type i ) {

        if(indices.at(i) >= _data.size() - (indices.size() - i)) {
            if( i!=0 && nextlevel(i-1) ) {
                indices.at(i) = indices.at(i-1) +1;
                return true;
            } else
                return false;
        } else {
            indices.at(i)++;
            return true;
        }
    }

public:
    typedef T value_type;

    /**
     * @brief KofNvector
     * @param _data The std::vector to draw from
     * @param k number of elemets to draw each time
     */
    NchooseK( const std::vector<T>& data, index_type k): _data(data), _done(false) {

        if(k>_data.size()) {
            k=0;
            _done=true;
        }

        indices.resize(k);

        for(index_type i=0;i<k; ++i) {
            indices.at(i) = i;
        }
    }

    /// \todo make a move constructor

    /**
     * @brief Access the ith element of the currently drawn combination
     * @param i
     * @return A reference to the element
     * @note read-only. If makes no sense ot modify elements in a combination.
     */
    const T& at( const size_t i ) const { return _data.at(indices.at(i)); }

    /**
     * @brief Generate the next combination.
     * @return false if no more combinations to do.
     */
    bool next() {
        if(k()==0) {
            _done = true;
            return false;
        }
        bool res = nextlevel(indices.size()-1);
        _done = !res;
        return res;
    }


    /**
     * @brief Iterator over the elements of a combination. Used to iterate over the current combination.
     */
    class const_iterator : public std::iterator<std::forward_iterator_tag, T>
    {
        index_list::const_iterator index;
        const NchooseK<T>& v;
    public:

        const_iterator( const NchooseK<T>& c, const index_list::const_iterator& i) : index(i), v(c) {}

        const_iterator(const const_iterator& mit) : index(mit.index), v(mit.v) {}

        const_iterator& operator++() {++index;return *this;}

        bool operator==(const const_iterator& rhs) {return index==rhs.index;}

        bool operator!=(const const_iterator& rhs) {return index!=rhs.index;}

        const T& operator*() const { return v._data.at(*index); }

        typedef T value_type;
    };

    const_iterator begin() const { return !_done ? const_iterator(*this, indices.begin()): end(); }
    const_iterator end() const { return const_iterator(*this, indices.end()); }

    const index_list& Indices() const { return indices; }

    NchooseK<T>& operator++ () { next(); return *this;}

    bool done() const { return _done; }

    std::size_t size() const {
        // it's n choose k, binomial coeefficient = fac(n)/fac(
        // use gamma as factorial function
        const auto fac = [] (std::size_t i) { return std::tgamma(i+1); };
        // note: this might fail for larger numbers, as factorials get big quite quickly
        // would be better to calculate expansion n(n-1)(n-2)...(n-k+1)/(k(k-1)(k-2)...1)
        return fac(n())/(fac(k())*fac(n()-k()));
    }

    /**
     * @brief k
     * @return number of elements to choose
     */
    std::size_t k() const { return indices.size(); }

    /**
     * @brief n
     * @return number of elements to choose from
     */
    std::size_t n() const { return _data.size();}
};

template <typename T>
NchooseK<T> makeCombination( const std::vector<T>& data, const unsigned int k) {
    return NchooseK<T>(data,k);
}

}}} // namespace ant::analysis::utils
