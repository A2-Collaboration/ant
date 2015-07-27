#pragma once

#include <vector>
#include <memory>
#include <type_traits>

#include "base/std_ext.h"

namespace ant {
namespace std_ext {

/**
 * @brief The mapped_items struct
 */


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

    static_assert(std::is_integral<Key>::value, "Key for mapped_vectors must be integral type");

    void init(Key maxKey) {
        storage.resize(maxKey+1);
        clear();
    }

    void clear() {
        for(auto key : keys)
            storage[key]->second.clear();
        keys.clear();
    }

    void add_item(const Key& key, const Value& value) {
        // return silently if we encounter
        // a key which is too large for the storage
        if(key>=storage.size())
            return;

        // retrieve an already existing key
        auto& ptr = storage[key];
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
            return *(storage_ptr->at(*it_key));
        }
    private:
        friend class mapped_vectors;

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
