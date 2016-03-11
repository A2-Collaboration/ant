#pragma once

#include <vector>
#include <memory>
#include <type_traits>

#include "base/std_ext/memory.h" // make_unique

namespace ant {
namespace std_ext {

template<class E, class Enable = void>
struct to_integral_helper
{
    static E inner(E e)
    {
        return e;
    }
};

template<typename E>
struct to_integral_helper<E, typename std::enable_if<std::is_enum<E>::value>::type>
{
    static typename std::underlying_type<E>::type inner(E e)
    {
        return static_cast<typename std::underlying_type<E>::type>(e);
    }
};

template<typename E>
auto to_integral(E e) -> decltype(to_integral_helper<E>::inner(e))
{
    return to_integral_helper<E>::inner(e);
}


/**
 * @brief The mapped_vectors struct provides fast access to a collection of vectors
 */
template<typename Key, typename Value>
struct mapped_vectors {

private:
    using storage_item_t = std::pair<Key, std::vector<Value> >;
    using storage_item_ptr_t = std::unique_ptr<storage_item_t>;
    using storage_t = std::vector< storage_item_ptr_t >;
    storage_t storage;
    using keys_t = std::vector<Key>;
    keys_t keys;

public:

    void clear() {
        for(auto key : keys) {
            const auto key_u = to_integral(key);
            storage[key_u]->second.resize(0);
        }
        keys.resize(0);
    }

    std::size_t size() const {
        return keys.size();
    }

    void add_item(const Key& key, const Value& value) {
        // return silently if we encounter
        // a key which is too large for the storage
        const auto key_u = to_integral(key);
        if(key_u>=storage.size()) {
            storage.resize(key_u+1);
        }

        // retrieve an already existing key
        auto& ptr = storage[key_u];
        if(ptr==nullptr) {
            ptr = std_ext::make_unique< std::pair<Key, std::vector<Value>> >(make_pair(key, std::vector<Value>()));
        }

        std::vector<Value>& values = ptr->second;

        // first time something is added, remember the key
        // for faster iteration later
        if(values.empty())
            keys.push_back(key);

        // add the value finally
        ptr->second.push_back(value);
    }

    std::vector<Value> get_item(const Key& key) const {
        // return silently if we encounter
        // a key which is too large for the storage
        const auto key_u = to_integral(key);
        if(key_u>=storage.size())
            return {};
        auto& ptr = storage[key_u];
        if(ptr==nullptr)
            return {};
        return ptr->second;
    }

    friend class const_iterator;

    class const_iterator : public std::iterator<std::forward_iterator_tag, storage_item_t> {
    public:
        using typename std::iterator<std::forward_iterator_tag, storage_item_t>::reference;
        using typename std::iterator<std::forward_iterator_tag, storage_item_t>::pointer;

        using const_reference =  typename std::add_const<reference>::type;
        using const_pointer =  typename std::add_const<pointer>::type;

        bool operator==(const const_iterator& rhs) const {
            return it_key == rhs.it_key;
        }
        bool operator!=(const const_iterator& rhs) const {
            return !(*this == rhs);
        }

        const_iterator& operator++() {
            ++it_key; return *this;
        }

        const_reference operator*() const {
            return *(storage_ptr->at(to_integral(*it_key)));
        }

        const_pointer operator->() const {
            return storage_ptr->at(to_integral(*it_key)).get();
        }

    private:
        friend struct mapped_vectors;

        using storage_ptr_t = typename mapped_vectors::storage_t;
        using it_key_t = typename keys_t::const_iterator;

        const_iterator(it_key_t it_key_, const storage_ptr_t* storage_ptr_) :
            storage_ptr(storage_ptr_),
            it_key(it_key_)
        {}

        const storage_ptr_t* storage_ptr;
        it_key_t it_key;
    };

    const_iterator begin() const {
        return const_iterator(keys.cbegin(), std::addressof(storage));
    }

    const_iterator end() const {
        return const_iterator(keys.cend(), std::addressof(storage));
    }
};

}}
