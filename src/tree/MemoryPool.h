#pragma once

#include "base/std_ext.h" // for make_unique

#include <memory>
#include <forward_list>


namespace ant {

template<class T>
struct MemoryPool {

    class Item {
        friend class MemoryPool;
        std::unique_ptr<T> ptr;
        Item(std::unique_ptr<T> item_) :
            ptr(move(item_)) {}
    public:
        T&  operator*() const { return *ptr.get(); }
        T* operator->() const noexcept { return ptr.get(); }
        ~Item() {
            if(ptr == nullptr)
               return;
            std::cout << "Return to pool" << std::endl;
            MemoryPool::ReturnToPool(move(ptr));
        }
        Item(Item&&) = default;
        Item& operator=(Item&&) = default;
    };

    static Item Get() {
        if(items.empty()) {
            std::cout << "Created new pool item" << std::endl;
            return Item(std_ext::make_unique<T>());
        }
        else {
            items.front()->Clear();
            Item item(move(items.front()));
            items.pop_front();
            std::cout << "REUSED pool item" << std::endl;
            return item;
        }
    }

    MemoryPool() = delete;
    MemoryPool(const MemoryPool&) = delete;
    MemoryPool& operator=(const MemoryPool&) = delete;
    MemoryPool(MemoryPool&&) = delete;
    MemoryPool& operator=(MemoryPool&&) = delete;

private:
    static std::forward_list<std::unique_ptr<T>> items;
    static void ReturnToPool(std::unique_ptr<T> ptr) {
        items.push_front(move(ptr));
    }
};

}
