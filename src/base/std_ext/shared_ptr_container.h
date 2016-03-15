#pragma once

#include <vector>
#include <memory>

namespace ant {
namespace std_ext {

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

    shared_ptr_iterator_t() {}

    // convert from non-const to const iterator
    template<class other_it_t>
    shared_ptr_iterator_t(const shared_ptr_iterator_t<other_it_t, false>& i) : it(i.it) {}

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
    std::shared_ptr<const value_type> get_const() const {
        return *it;
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
    using c_t = Container<std::shared_ptr<T>, std::allocator<std::shared_ptr<T>>>;
    c_t c;

public:

    using const_iterator = shared_ptr_iterator_t<typename c_t::const_iterator>;
    using iterator = shared_ptr_iterator_t<typename c_t::iterator>;
    using items_t = std::initializer_list<const_iterator>;

    shared_ptr_container() = default;

    template<class it_t>
    shared_ptr_container(std::initializer_list< shared_ptr_iterator_t<it_t> > items)
    {
        // need to convert iterators from shared_ptr here
        // since only this class has access to shared_ptr iterator
        for(const auto& i : items)
            c.emplace_back(*i.it);
    }

    template<class it_t>
    shared_ptr_container(shared_ptr_iterator_t<it_t> first, shared_ptr_iterator_t<it_t> last)
    {
        // need to convert iterators from shared_ptr here
        // since only this class has access to shared_ptr iterator
        for(auto i = first; i != last; ++i)
            c.emplace_back(*i.it);
    }


    // copy/move
    shared_ptr_container(const shared_ptr_container&) = delete;
    shared_ptr_container& operator=(const shared_ptr_container&) = delete;
    shared_ptr_container(shared_ptr_container&&) = default;
    shared_ptr_container& operator=(shared_ptr_container&&) = default;

    const_iterator begin() const noexcept
    {
        return c.begin();
    }

    const_iterator end() const noexcept
    {
        return c.end();
    }

    iterator begin() noexcept
    {
        return c.begin();
    }

    iterator end() noexcept
    {
        return c.end();
    }

    template<class... Args>
    void emplace_back(Args&&... args) {
        c.emplace_back(std::make_shared<T>(std::forward<Args>(args)...));
    }

    void emplace_back(const_iterator item) {
        c.emplace_back(*item.it);
    }

    template<class it_t>
    iterator insert(const_iterator pos, it_t first, it_t last) {
        return c.insert(pos.it, first.it, last.it);
    }


    const T& operator[] (typename c_t::size_type i) const noexcept
    {
        return *c[i];
    }

    const T& at(typename c_t::size_type i) const {
        return *c.at(i);
    }

    const T& back() const {
        return *c.back();
    }

    T& back() {
        return *c.back();
    }

    const T& front() const {
        return *c.front();
    }

    T& front() {
        return *c.front();
    }


    iterator erase(const iterator& it) noexcept
    {
        return c.erase(it.it);
    }

    typename c_t::size_type size() const noexcept {
        return c.size();
    }

    bool empty() const noexcept {
        return c.empty();
    }

    void clear() noexcept {
        c.clear();
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(c);
    }
};


}} // namespace ant::std_ext