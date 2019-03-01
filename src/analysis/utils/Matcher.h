#pragma once

#include <list>
#include <vector>
#include <algorithm>

#include "base/interval.h"
#include "base/std_ext/math.h"
#include <limits>



namespace ant {
namespace analysis {
namespace utils {

/**
 * @gbrief Matcher Function Template: Match by simple distance |a-b|
 * @param a
 * @param b
 * Use for scalar types
 */
template <typename T>
static double matchDistance( const T& a, const T& b ) {
    return fabs( double(a - b) );
}

template <typename T1, typename T2>
struct scored_match {
    double score;
    T1 a;
    T2 b;

    bool operator< (const scored_match<T1,T2>& rhs ) const { return score < rhs.score; }

    std::pair<T1,T2> getPair() const { return std::pair<T1,T2>(a,b); }
};

/**
 * @brief Particle Matcher
 * @param list1 vector of elements 1
 * @param list2 vector of elements 2
 * @param f
 * @param score_window
 *
 * Matches the elements in list1 and list2 by a score calculated by a MatchFunction.
 * returns a list of best-matching pairs (minimal score),
 * every element of list1,list2 occurs only once in the pair list
 *
 * The Matcher Function can be a lambda function
 */
template <class MatchFunction, typename List1, typename List2>
std::list< scored_match<typename List1::value_type, typename List2::value_type> >
    match1to1( const List1& list1,
               const List2& list2,
               MatchFunction f,
               const ant::IntervalD& score_window=ant::IntervalD(-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()) )
{

    typedef typename List1::value_type T1;
    typedef typename List2::value_type T2;
    typedef std::list< scored_match<T1, T2> > scorelist;


    scorelist scores;

    // build pairs and calculate scores
    for( const auto& i : list1 )
        for( const auto& j : list2 ) {
            const double score = f(i,j);
            if( score_window.Contains(score)) {
                scored_match<T1,T2> s = {score,i,j};
                scores.emplace_back(s);
            }
        }

    scores.sort();

    // find best matching pairs
    auto i = scores.begin();
    while(  i != scores.end() ) {

        auto j = i;
        ++j;

        //remove all matches that include one of the two just matched particles
        while( j!=scores.end() ) {

            if( j->a == i->a || j->b == i->b)
                j = scores.erase(j);
            else
                ++j;
        }

        ++i;
    }

    return scores;
}

template <typename T1, typename T2>
T2 FindMatched(const std::list<utils::scored_match<T1,T2>>& l, const T1& f) {
    for( const auto& i : l ) {
        if( i.a == f) {
            return i.b;
        }
    }
    return T2();
}

// this method returns a list of particles
// which are not in the in the scored list
template <typename T1, typename T2>
std::vector<T2> FindUnmatched(const std::list<utils::scored_match<T1,T2>>& l, const std::vector<T2>& f) {
    std::vector<T2> unmatched;
    for( const auto& i : f ) { // loop over std::vector<T2>
        bool found = false;
        for( const auto& j : l ) { // loop over list of scored_match<T1,T2>
            if( j.b == i) {
                found = true;
                continue;
            }
        }
        if (found == false) {
            unmatched.push_back(i);
        }
    }
    return unmatched;
};



struct matchpair {
    size_t a;
    size_t b;
    double dist;
    bool matched = false;

    matchpair(size_t index_a, size_t index_b, double distance) noexcept:
        a(index_a),
        b(index_b),
        dist(distance) {}

    bool operator> (const matchpair& o) const noexcept {
        return dist > o.dist;
    }

    bool operator< (const matchpair& o) const noexcept {
        return dist < o.dist;
    }
};

template <class MatchFunction, typename List1, typename List2>
std::vector<matchpair>match2(const List1& list1,
               const List2& list2,
               MatchFunction f)
{
    std::list<matchpair> pairs;

    for(size_t i=0; i<list1.size(); ++i) {
        for(size_t j=0; j<list2.size(); ++j) {
            pairs.emplace_back(matchpair(i,j, f(list1.at(i), list2.at(j))));
        }
    }

    std::vector<matchpair> bestlist;
    bestlist.reserve(list1.size());
    for(size_t i=0;i<list1.size(); ++i) {
        bestlist.emplace_back(matchpair(i, std::numeric_limits<size_t>::max(), std_ext::inf));
    }

    // save best match for every element in A,
    // regardless of other matches.
    for(const auto& m : pairs) {
        auto& b = bestlist.at(m.a);
        if(m.dist < b.dist) {
            b = m;
        }
    }

    pairs.sort();

    {

        while(!pairs.empty()) {
            auto i = pairs.begin();

            auto& b = bestlist.at(i->a);
            b = *i;

            b.matched = true;

            while(i!=pairs.end()) {
                if(i->a == b.a || i->b == b.b) {
                    i = pairs.erase(i);
                } else {
                    ++i;
                }
            }
        }
    }

    return bestlist;
}




}}} // namespace ant::analysis::utils
