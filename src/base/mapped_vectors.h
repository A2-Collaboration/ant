#pragma once

#include <vector>
#include <memory>
#include <type_traits>

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
    using storage_item_t = std::unique_ptr< std::pair<Key, std::vector<Value> > >;
    using storage_t = std::vector< storage_item_t >;
    storage_t storage;
    using keys_t = std::vector<Key>;
    keys_t keys;

public:
//    using is_enum = std::is_enum<Key>;
//    using Key_enum = typename std::underlying_type<Key>::type;
//    using Key_u = typename std::conditional< is_enum::value, typename std::underlying_type<Key>::type, Key >::type;

//    static_assert(std::is_integral< Key_u >::value, "Key for mapped_vectors must be integral type");

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
            ptr = make_unique< std::pair<Key, std::vector<Value>> >(make_pair(key, std::vector<Value>()));
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

    using A = std::allocator< std::pair<Key, std::vector<Value>> >;

    class const_iterator {
    public:
        typedef typename A::difference_type difference_type;
        typedef typename A::value_type value_type;
        typedef typename A::reference const_reference;
        typedef typename A::pointer const_pointer;
        typedef std::forward_iterator_tag iterator_category;

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
