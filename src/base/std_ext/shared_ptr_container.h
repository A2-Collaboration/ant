#pragma once

#include <vector>
#include <memory>

namespace ant {
namespace std_ext {

template<typename T, class it_t>
struct shared_ptr_iterator_t : std::iterator<std::bidirectional_iterator_tag, T>
{
    static_assert(std::is_same<typename it_t::value_type, std::shared_ptr<T>>::value,
                  "shared_ptr_iterator_t must be used with iterator pointing to shared_ptr<T>");

    using base_type = typename std::iterator<std::bidirectional_iterator_tag, T>;
    using typename base_type::reference;
    using typename base_type::pointer;
    using const_reference =  typename std::add_const<reference>::type;
    using const_pointer =  typename std::add_const<pointer>::type;

    bool operator==(const shared_ptr_iterator_t& rhs) const {
        return it == rhs.it;
    }
    bool operator!=(const shared_ptr_iterator_t& rhs) const {
        return !(*this == rhs);
    }
    shared_ptr_iterator_t& operator++() {
        ++it; return *this;
    }
    shared_ptr_iterator_t& operator--() {
        --it; return *this;
    }
    const_reference operator*() const {
        return **it;
    }
    const_pointer operator->() const {
        return std::addressof(*this);
    }
    std::shared_ptr<const T> get_const() const {
        return *it;
    }

private:
    template<typename, template<class, class> class>
    friend class shared_ptr_container;

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

    using const_iterator = shared_ptr_iterator_t<T, typename c_t::const_iterator>;

    shared_ptr_container() = default;

    template<class it_t>
    shared_ptr_container(const std::initializer_list< shared_ptr_iterator_t<T, it_t> >& items)
    {
        // need to convert iterators from shared_ptr here
        // since only this class has access to shared_ptr iterator
        for(const auto& i : items)
            c.emplace_back(*i.it);
    }

    // copy/move
    shared_ptr_container(const shared_ptr_container&) = delete;
    shared_ptr_container& operator=(const shared_ptr_container&) = delete;
    shared_ptr_container(shared_ptr_container&&) = delete;
    shared_ptr_container& operator=(shared_ptr_container&&) = delete;

    const_iterator begin() const
    {
        return c.begin();
    }

    const_iterator end() const
    {
        return c.end();
    }

    template<typename... Args>
    void emplace_back(Args&&... args)
    {
        c.emplace_back(std::make_shared<T>(std::forward<Args>(args)...));
    }

    const T& operator[] (typename const_iterator::difference_type i) const
    {
        return **std::next(c.begin(), i);
    }

    const_iterator erase(const const_iterator& it)
    {
        // need explicit cast to const_iterator since erase returns normal iterator
        return const_iterator(c.erase(it.it));
    }
};


}} // namespace ant::std_ext