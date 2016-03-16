#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <algorithm>

namespace ant {
namespace std_ext {

template<class T>
struct cc_shared_ptr {

    using element_type = const T;

    cc_shared_ptr() = default;
    cc_shared_ptr(const cc_shared_ptr&) = default;
    cc_shared_ptr& operator=(const cc_shared_ptr&) = default;
    cc_shared_ptr(cc_shared_ptr&&) = default;
    cc_shared_ptr& operator=(cc_shared_ptr&&) = default;

    cc_shared_ptr(std::nullptr_t) : cc_shared_ptr() {}
    cc_shared_ptr(const std::shared_ptr<T>& ptr_) : ptr(ptr_) {}

    const T& operator*() const {
        return *ptr;
    }
    const T* operator->() const {
        return ptr.operator->();
    }
    bool operator==(const cc_shared_ptr& other) const {
        return ptr == other.ptr;
    }
    explicit operator bool() const {
        return static_cast<bool>(ptr);
    }
    // make it cerealizable
    template<class Archive>
    void serialize(Archive& archive) { archive(ptr); }
private:
    std::shared_ptr<T> ptr;
};

template<class it_t, class T = decltype(*std::declval<it_t>())>
constexpr bool
is_const_iterator() {
    return  ! std::is_assignable <
            decltype( *std::declval<it_t>() ),
            T
            >::value;
}


template<class it_t, bool is_const = is_const_iterator<it_t>()>
struct shared_ptr_iterator_t
{
    // assume that it_t is iterator over shared_ptr, which defines element_type
    using value_type = typename it_t::value_type::element_type;
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = typename it_t::difference_type;
    using pointer =  typename std::conditional<is_const, const value_type*, value_type*>::type;
    using reference =  typename std::conditional<is_const, const value_type&, value_type&>::type;

    // convert from non-const to const iterator
    template<typename other_it_t>
    shared_ptr_iterator_t(const shared_ptr_iterator_t<other_it_t, false>& i) :
        shared_ptr_iterator_t(i.it) {}

    shared_ptr_iterator_t& operator++() {
        ++it; return *this;
    }
    shared_ptr_iterator_t& operator--() {
        --it; return *this;
    }
    reference operator*() const {
        return **it;
    }
    pointer operator->() const {
        return std::addressof(**it);
    }
    cc_shared_ptr<value_type> get_const() const {
        return *it;
    }

    operator cc_shared_ptr<value_type>() const {
        return get_const();
    }

    friend bool operator==(const shared_ptr_iterator_t& x,
                           const shared_ptr_iterator_t& y)
    {
        return x.it == y.it;
    }

    friend bool operator!=(const shared_ptr_iterator_t& x,
                           const shared_ptr_iterator_t& y)
    {
        return x.it != y.it;
    }

private:
    template<typename, template<class, class> class>
    friend class shared_ptr_container;

    template<class, bool>
    friend class shared_ptr_iterator_t;

    shared_ptr_iterator_t(const it_t& it_) : it(it_) {}
    it_t it;
};

template<typename T, template<class, class> class Container = std::vector>
class shared_ptr_container
{
private:
    // shortcut to use container with std::allocator
    template<typename t, template<class, class> class container>
    using container_t = container<t, std::allocator<t>>;

    using c_t = container_t<std::shared_ptr<T>, Container>;
    c_t c;

public:
    using size_type = typename c_t::size_type;
    using const_iterator = shared_ptr_iterator_t<typename c_t::const_iterator>;
    using iterator = shared_ptr_iterator_t<typename c_t::iterator>;

    template<class it_t>
    using non_const_iterator = shared_ptr_iterator_t<it_t, false>;

    shared_ptr_container() = default;

    template<class it_t>
    shared_ptr_container(std::initializer_list<non_const_iterator<it_t>> items)
    {
        for(auto& item: items)
            c.emplace_back(*item.it);
    }

    template<class it_t>
    shared_ptr_container(non_const_iterator<it_t> first, non_const_iterator<it_t> last) :
        c(first.it, last.it)
    {}

    // make it cerealizable
    template<class Archive>
    void serialize(Archive& archive) { archive(c); }

    // make it movable only
    shared_ptr_container(const shared_ptr_container&) = delete;
    shared_ptr_container& operator=(const shared_ptr_container&) = delete;
    shared_ptr_container(shared_ptr_container&&) = default;
    shared_ptr_container& operator=(shared_ptr_container&&) = default;

    // access methods
    const_iterator begin() const noexcept { return c.begin(); }
          iterator begin()       noexcept { return c.begin(); }

    const_iterator end() const noexcept { return c.end(); }
          iterator end()       noexcept { return c.end(); }

    const T& back() const { return *c.back(); }
          T& back()       { return *c.back(); }

    const T& front() const { return *c.front(); }
          T& front()       { return *c.front(); }

    const T& operator[] (typename c_t::size_type i) const noexcept
    {
        return *c[i];
    }

    const T& at(typename c_t::size_type i) const
    {
        return *c.at(i);
    }

    cc_shared_ptr<T> get_const(size_type i) const noexcept
    {
        return c[i];
    }

    template<template<class, class> class container = Container>
    container_t<cc_shared_ptr<T>, container> get_const_list(
            std::function<bool(const T&)> filter
            = [] (const T&) {return true;}) const {
        container_t<cc_shared_ptr<T>, container> tmp_c;
        for(auto& item : c) {
            if(filter(*item))
                tmp_c.emplace_back(item);
        }
        return tmp_c;
    }

    size_type size() const noexcept { return c.size(); }

    bool empty() const noexcept { return c.empty(); }

    void clear() noexcept { c.clear(); }

    template<class... Args>
    void emplace_back(Args&&... args)
    {
        c.emplace_back(std::make_shared<T>(std::forward<Args>(args)...));
    }

    template<class it_t>
    void push_back(non_const_iterator<it_t> item)
    {
        c.emplace_back(*item.it);
    }

    template<class it_t>
    void insert(iterator pos,
                non_const_iterator<it_t> first,
                non_const_iterator<it_t> last)
    {
        c.insert(pos.it, first.it, last.it);
    }

    iterator erase(const iterator& it) noexcept
    {
        return c.erase(it.it);
    }
};


}} // namespace ant::std_ext