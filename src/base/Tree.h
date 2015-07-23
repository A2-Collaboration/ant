#pragma once

#include <memory>
#include <list>

namespace ant {

template <typename T>
class Tree {
protected:
    T data;

    using snode_t    = std::shared_ptr<Tree>;
    using wnode_t    = std::weak_ptr<Tree>;
    using snodelist_t = std::list<snode_t>;

    wnode_t partent;
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
        return std::move(n);
    }

    bool isRoot() const { return partent.expired(); }
    bool isLeaf() const { return daughters.empty(); }

    T& operator*() { return data; }
    const T& operator *() const { return data; }

    const T& get() const { return data; }

    snode_t GetParent() { return partent.lock(); }
    snode_t Self() { return self.lock(); }
    const snodelist_t& Daughters() const { return daughters; }

    void Unlink() {
        if(auto p = GetParent()) {
            p->remove_daughter(this);
            partent.reset();
        }
    }

    template <typename ... args_t>
    snode_t createDaughter(args_t&&... args) {
        auto n = makeNode(std::forward<args_t>(args)...);
        n->partent = Self();
        daughters.emplace_back(n);
        return std::move(n);
    }

    void setParent(snode_t& p) {
        if(!isRoot())
            Unlink();
        partent = p;
        p->daughters.emplace_back(Self);
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
        maplevel(
                    [&d] (const T&, const size_t& level) {
            if(level > d) {d=level;}}
    );
        return d +1;
    }

};

}
