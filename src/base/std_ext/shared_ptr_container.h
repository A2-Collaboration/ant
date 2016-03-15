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

    using const_reference =  const T&;
    using const_pointer = const T*;

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
        return std::addressof(**it);
    }
    std::shared_ptr<const T> get_const() const {
        return *it;
    }

    explicit operator bool() const {
        return it != it_end;
    }

    shared_ptr_iterator_t() {}

private:
    template<typename, template<class, class> class>
    friend class shared_ptr_container;

    shared_ptr_iterator_t(const it_t& it_, const it_t& it_end_) : it(it_), it_end(it_end_) {}
    it_t it;
    it_t it_end;

};

template<typename T, template<class, class> class Container = std::vector>
class shared_ptr_container
{
private:
    using c_t = Container<std::shared_ptr<T>, std::allocator<std::shared_ptr<T>>>;
    c_t c;

public:

    using const_iterator = shared_ptr_iterator_t<T, typename c_t::const_iterator>;
    using items_t = std::initializer_list<const_iterator>;

    shared_ptr_container() = default;

    template<class it_t>
    shared_ptr_container(std::initializer_list< shared_ptr_iterator_t<T, it_t> > items)
    {
        // need to convert iterators from shared_ptr here
        // since only this class has access to shared_ptr iterator
        for(const auto& i : items)
            c.emplace_back(*i.it);
    }

    template<class it_t>
    shared_ptr_container(shared_ptr_iterator_t<T, it_t> first, shared_ptr_iterator_t<T, it_t> last)
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
        return const_iterator(c.begin(), c.end());
    }

    const_iterator end() const noexcept
    {
        return const_iterator(c.end(), c.end());
    }

    template<class... Args>
    void emplace_back(Args&&... args) {
        c.emplace_back(std::make_shared<T>(std::forward<Args>(args)...));
    }

    void emplace_back(const_iterator item) {
        c.emplace_back(*item.it);
    }

    const_iterator insert(const_iterator pos, const_iterator first, const_iterator last) {
        auto it = c.insert(pos.it, first.it, last.it);
        return const_iterator(it, c.end());
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

    const_iterator erase(const const_iterator& it) noexcept
    {
        // need explicit cast to const_iterator since erase returns normal iterator
        auto it_ = c.erase(it.it);
        return const_iterator(it_, c.end());
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