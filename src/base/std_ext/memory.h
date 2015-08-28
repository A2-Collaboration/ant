#pragma once

#include <memory>
#include <list>

namespace ant {
namespace std_ext {

template<typename T, typename Base>
void AddToSharedPtrList(const std::shared_ptr<Base> base, std::list< std::shared_ptr<T> >& list) {
    const auto& ptr = std::dynamic_pointer_cast<T, Base>(base);
    if(ptr != nullptr)
        list.emplace_back(std::move(ptr));
}

template<typename From, typename To>
inline std::unique_ptr<To> static_cast_uptr(std::unique_ptr<From>&& ptr) {
  static_assert(std::is_base_of<To, From>::value, "Type From does not derive from type To");
  return std::unique_ptr<To>(static_cast<To*>(ptr.release()));
}

// copied from C++14 draft
// https://isocpp.org/files/papers/N3656.txt
// once we use C++14 by default, we can remove this here
template<class T> struct _Unique_if {
  typedef std::unique_ptr<T> _Single_object;
};

template<class T> struct _Unique_if<T[]> {
  typedef std::unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N> struct _Unique_if<T[N]> {
  typedef void _Known_bound;
};

template<class T, class... Args>
typename _Unique_if<T>::_Single_object
make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class T>
typename _Unique_if<T>::_Unknown_bound
make_unique(size_t n) {
  typedef typename std::remove_extent<T>::type U;
  return std::unique_ptr<T>(new U[n]());
}

template<class T, class... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args&&...) = delete;

}} // namespace ant::std_ext
