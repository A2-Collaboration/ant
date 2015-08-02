#pragma once

#include "base/std_ext.h" // for make_unique

#include <memory>
#include <forward_list>


namespace ant {

template<class T>
struct MemoryPool {

    class Item {
        friend class MemoryPool;
        MemoryPool* const pool;
        std::unique_ptr<T> ptr;
        Item(MemoryPool* pool_, std::unique_ptr<T> ptr_) :
            pool(pool_),
            ptr(std::move(ptr_))
        {}
    public:
        T* get() const { return ptr.get(); }
        T&  operator*() const { return *get(); }
        T* operator->() const noexcept { return get(); }
        ~Item() {
            if(ptr == nullptr || pool == nullptr)
               return;
            pool->ReturnToPool(std::move(ptr));
        }
        Item(Item&&) = default;
        Item& operator=(Item&&) = default;
    };

    static Item Get() {
        static MemoryPool m;
        if(m.items.empty()) {
            return Item(std::addressof(m), std_ext::make_unique<T>());
        }
        else {
            m.items.front()->Clear();
            Item item(std::addressof(m), std::move(m.items.front()));
            m.items.pop_front();
            return item;
        }
    }

    MemoryPool() = default;
    MemoryPool(const MemoryPool&) = delete;
    MemoryPool& operator=(const MemoryPool&) = delete;
    MemoryPool(MemoryPool&&) = delete;
    MemoryPool& operator=(MemoryPool&&) = delete;

private:
    std::forward_list<std::unique_ptr<T>> items;
    void ReturnToPool(std::unique_ptr<T> ptr) {
        items.push_front(std::move(ptr));
    }
};

}
