#pragma once

#include <memory>
#include <list>
#include <algorithm>

namespace ant {

template <typename T>
class Tree {
protected:
    T data;

    using snode_t    = std::shared_ptr<Tree>;
    using wnode_t    = std::weak_ptr<Tree>;
    using snodelist_t = std::list<snode_t>;

    wnode_t parent;
    wnode_t self;
    snodelist_t daughters;

    void remove_daughter(Tree* t) {
        daughters.remove_if( [t] (snode_t& n) { return n.get() == t;});
    }


    Tree() {}
    Tree(const T& data_): data(data_) {}

public:

    template <typename ... args_t>
    static snode_t makeNode(args_t&&... args) {
        auto n = std::shared_ptr<Tree>(new Tree(std::forward<args_t>(args)...));
        n->self = n;
        return n;
    }

    bool isRoot() const { return parent.expired(); }
    bool isLeaf() const { return daughters.empty(); }

    T& operator*() { return data; }
    const T& operator *() const { return data; }

    const T& get() const { return data; }

    snode_t GetParent() { return parent.lock(); }
    snode_t Self() { return self.lock(); }
    const snodelist_t& Daughters() const { return daughters; }

    void Unlink() {
        if(auto p = GetParent()) {
            p->remove_daughter(this);
            parent.reset();
        }
    }

    template <typename ... args_t>
    snode_t createDaughter(args_t&&... args) {
        auto n = makeNode(std::forward<args_t>(args)...);
        n->parent = Self();
        daughters.emplace_back(n);
        return n;
    }

    void setParent(snode_t& p) {
        if(!isRoot())
            Unlink();
        parent = p;
        p->daughters.emplace_back(Self());
    }

    template <typename F>
    void map(F function) const {
        function(get());
        for(auto& node: daughters) {
            node->map(function);
        }
    }

    template <typename F>
    void maplevel(F function, size_t level=0) const {
        function(get(), level);
        for(auto& node: daughters) {
            node->maplevel(function, level+1);
        }
    }

    size_t size() const {
        size_t s = 0;
        map([&s] (const T&) { s++;});
        return s;
    }

    size_t depth() const {
        size_t d = 0;
        maplevel([&d] (const T&, const size_t& level) {
            if(level>d) d=level;
        }
        );
        return d +1;
    }

    template<typename Compare>
    void sort(Compare comp) {
        // sort daughters first (depth-first recursion)
        std::for_each(daughters.begin(), daughters.end(),
                      [comp] (snode_t d) { d->sort(comp); });

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
    }

};

}
