#ifndef TREE_H
#define TREE_H

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
    snodelist_t daughters;

    void remove_daughter(Tree* t) {
        daughters.remove_if( [t] (snode_t& n) { return n.get() == t;});
    }

public:

    Tree() {}
    Tree(const T& data_): data(data_) {}

    static snode_t makeNode(const T& data) { return std::move(std::make_shared<Tree>(data)); }

    static snode_t makeNode(snode_t& p, const T& data) {
        auto n =  makeNode(data);
        link(p,n);
        return std::move(n);
    }


    bool isRoot() const { return partent.expired(); }
    bool isLeaf() const { return daughters.empty(); }

    T& operator*() { return data; }
    const T& operator *() const { return data; }

    const T& get() const { return data; }

    static void link( snode_t& p, snode_t& daughter) {

        if(!p || !daughter)
            return;

        if( !daughter->isRoot() ) {
            daughter->Unlink();
        }
        daughter->partent = p;
        p->daughters.emplace_back(daughter);
    }

    snode_t GetParent() { return partent.lock(); }

    void Unlink() {
        if(auto p = GetParent()) {
            p->remove_daughter(this);
            partent.reset();
        }
    }

    const snodelist_t& Daughters() const { return daughters; }

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

#endif
