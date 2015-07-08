#ifndef FORMAT_H
#define FORMAT_H

#include <vector>
#include "detail/format.h"

//// this code makes it possible to format a vector of doubles
//namespace fmt {
//namespace ant_detail {

//// pack of indices pattern with variadic templates
//template< std::size_t... Ns >
//struct indices {
//    typedef indices< Ns..., sizeof...( Ns ) > next;
//};

//template< std::size_t N >
//struct make_indices {
//    typedef typename make_indices< N - 1 >::type::next type;
//};

//template<>
//struct make_indices< 0 > {
//    typedef indices<> type;
//};

//template<typename T, size_t... Is>
//std::string format_vector_helper(StringRef format_str, const std::vector<T>& v, indices<Is...>) {
//  return format(format_str, v[Is]...);
//}


//} // namespace fmt::ant_detail

//template<typename T>
//std::string format_vector(StringRef format_str, const std::vector<T>& v) {
//  //return format(format_str, args);
//}

//} // namespace fmt

#endif // FORMAT_H
