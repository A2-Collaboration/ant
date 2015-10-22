#pragma once

#include <memory>
#include <list>
#include <algorithm>

namespace ant {

template<typename T>
class Tree {
    // make Tree of different types friends
    // needed for IsEqual()
    template<typename>
    friend class Tree;

protected:
    T data;

    using snode_t    = std::shared_ptr<Tree>;
    using wnode_t    = std::weak_ptr<Tree>;
    using snodelist_t = std::list<snode_t>;

    wnode_t parent;
    wnode_t self;
    snodelist_t daughters;
    bool is_sorted; // internal flag for IsEqual()

    void RemoveDaughter(Tree* t) {
        daughters.remove_if( [t] (snode_t& n) { return n.get() == t;});
    }


    Tree(T&& data_) : data(data_), is_sorted(false) {}

public:

    using type = T;

    template <typename ... args_t>
    static snode_t MakeNode(args_t&&... args) {
        // cannot use make_shared since protected constructor
        auto n = std::shared_ptr<Tree>(new Tree(T(std::forward<args_t>(args)...)));
        n->self = n;
        return n;
    }

    bool IsRoot() const { return parent.expired(); }
    bool IsLeaf() const { return daughters.empty(); }

    T& operator*() { return data; }
    const T& operator *() const { return data; }

    const T& Get() const { return data; }

    snode_t GetParent() const { return parent.lock(); }
    snode_t Self() const { return self.lock(); }
    const snodelist_t& Daughters() const { return daughters; }

    void Unlink() {
        if(auto p = GetParent()) {
            p->RemoveDaughter(this);
            p->is_sorted = false;
            parent.reset();
        }
    }

    template <typename ... args_t>
    snode_t CreateDaughter(args_t&&... args) {
        auto n = MakeNode(std::forward<args_t>(args)...);
        n->parent = Self();
        daughters.emplace_back(n);
        is_sorted = false;
        return n;
    }

    void SetParent(const snode_t& p) {
        if(!IsRoot())
            Unlink();
        parent = p;
        p->is_sorted = false;
        p->daughters.emplace_back(Self());
    }

    template <typename F>
    void Map(F function) const {
        function(Get());
        for(auto& node: daughters) {
            node->Map(function);
        }
    }

    template <typename F>
    void Map_level(F function, size_t level=0) const {
        function(Get(), level);
        for(auto& node: daughters) {
            node->Map_level(function, level+1);
        }
    }

    size_t Size() const {
        size_t s = 0;
        Map([&s] (const T&) { s++;});
        return s;
    }

    size_t Depth() const {
        size_t d = 0;
        Map_level([&d] (const T&, const size_t& level) {
            if(level>d) d=level;
        }
        );
        return d +1;
    }

    void Sort() {
        Sort(std::less<T>());
    }

    template<typename Compare>
    void Sort(Compare comp) {
        // sort daughters first (depth-first recursion)
        std::for_each(daughters.begin(), daughters.end(),
                      [comp] (snode_t d) { d->Sort(comp); });

        // then sort daughters itself
        daughters.sort([comp] (snode_t a, snode_t b) {
            auto a_less_b = comp(a->data, b->data);
            auto b_less_a = comp(b->data, a->data);
            auto equal = !a_less_b && !b_less_a;
            // sort by data first
            if(!equal)
                return a_less_b;
            // if equal, have a look at the sorted daughters itself
            if(a->daughters.size() != b->daughters.size())
                return a->daughters.size() < b->daughters.size();
            auto d_a = a->daughters.begin();
            auto d_b = b->daughters.begin();
            for(; d_a != a->daughters.end() && d_b != b->daughters.end(); ++d_a, ++d_b) {
                auto d_a_less_b = comp((*d_a)->data, (*d_b)->data);
                auto d_b_less_a = comp((*d_b)->data, (*d_a)->data);
                auto d_equal = !d_a_less_b && !d_b_less_a;
                if(!d_equal)
                    return d_a_less_b;
            }
            return false;
        });

        is_sorted = true;
    }

    template<typename U, typename Compare>
    bool IsEqual(const std::shared_ptr<Tree<U>>& other, Compare comp) const {
        // check daughters first (depth-first recursion)
        if(daughters.size() != other->daughters.size())
            return false;
        auto d = daughters.begin();
        auto d_other = other->daughters.begin();
        for(; d != daughters.end() && d_other != other->daughters.end(); ++d, ++d_other) {
            if(!(*d)->IsEqual(*d_other, comp))
                return false;
        }

        if(!is_sorted || !other->is_sorted)
            throw std::runtime_error("Can only compare sorted trees to each other");

        return comp(data, other->data);
    }

};

} // namespace ant
