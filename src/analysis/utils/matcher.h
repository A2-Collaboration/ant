#pragma once

#include <list>
#include <vector>
#include <algorithm>

#include "base/interval.h"
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

}}} // namespace ant::analysis::utils
