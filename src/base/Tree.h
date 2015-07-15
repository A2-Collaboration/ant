#ifndef TREE_H
#define TREE_H

#include <memory>
#include <list>

namespace ant {

template <typename T>
class Tree {
protected:
    T data;

    using node_t = std::shared_ptr<Tree>;
    using nodelist_t = std::list<node_t>;

    node_t partent;
    nodelist_t daughters;

    void remove_daughter(Tree* t) {
        daughters.remove_if( [t] (node_t& n) { return n.get() == t;});
    }

public:

    Tree() {}
    Tree(const T& data_): data(data_) {}

    static node_t makeNode(const T& data) { return std::make_shared<Tree>(data); }

    static node_t makeNode(node_t p, const T& data) {
        auto n =  std::make_shared<Tree>(data);
        link(p,n);
        return std::move(n);
    }


    bool isRoot() const { return partent.get() == nullptr; }
    bool isLeaf() const { return daughters.empty(); }

    T& operator*() { return data; }
    const T& operator *() const { return data; }

    const T& get() const { return data; }

    static void link( node_t& p, node_t& daughter) {
        if( !daughter->isRoot() ) {
            daughter->Unlink();
        }
        daughter->partent = p;
        p->daughters.push_back(daughter);
    }

    void Unlink() {
        if(!isRoot()) {
            partent->remove_daughter(this);
            partent.reset();
        }
    }

    const nodelist_t& Daughters() const { return daughters; }


};

}

#endif
