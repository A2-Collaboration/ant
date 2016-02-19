#pragma once

#include "base/cereal/access.hpp"

#include <memory>
#include <list>
#include <vector>
#include <algorithm>
#include <cassert>

namespace ant {

template<typename T>
class Tree {
    // make Tree of different types friends
    // needed for IsEqual()
    template<typename>
    friend class Tree;

    template<typename U>
    using snode_t    = std::shared_ptr<Tree<U>>;
public:

    using node_t = snode_t<T>;
    using type = T;

protected:
    T data;

    using wnode_t    = std::weak_ptr<Tree>;
    using snodelist_t = std::list<node_t>;

    wnode_t parent;
    wnode_t self;
    snodelist_t daughters;
    bool is_sorted; // internal flag for IsEqual()

    void RemoveDaughter(Tree* t) {
        daughters.remove_if( [t] (node_t& n) { return n.get() == t;});
    }


    Tree(T&& data_) : data(std::forward<T>(data_)), is_sorted(false) {}


    friend class cereal::access;
    Tree() {} // in order to construct shared_ptr of Tree<T>
    template <class Archive>
    void serialize(Archive& archive)
    {
        archive(data, parent, self, daughters, is_sorted);
    }

public:


    template <typename ... args_t>
    static node_t MakeNode(args_t&&... args) {
        // cannot use make_shared since protected constructor
        auto n = std::shared_ptr<Tree>(new Tree(T(std::forward<args_t>(args)...)));
        n->self = n;
        return n;
    }

    bool IsRoot() const { return parent.expired(); }
    bool IsLeaf() const { return daughters.empty(); }

    T& operator*() { return data; }
    const T& operator *() const { return data; }

    T& Get() { return data; }
    const T& Get() const { return data; }

    node_t GetParent() const { return parent.lock(); }
    node_t Self() const { return self.lock(); }
    const snodelist_t& Daughters() const { return daughters; }

    void Unlink() {
        if(auto p = GetParent()) {
            p->RemoveDaughter(this);
            p->is_sorted = false;
            parent.reset();
        }
    }

    template <typename ... args_t>
    node_t& CreateDaughter(args_t&&... args) {
        return AddDaughter(MakeNode(std::forward<args_t>(args)...));
    }

    node_t& AddDaughter(const node_t& n) {
        is_sorted = false;
        n->parent = Self();
        daughters.emplace_back(n);
        return daughters.back();
    }

    void SetParent(const node_t& p) {
        if(!IsRoot())
            Unlink();
        parent = p;
        p->is_sorted = false;
        p->daughters.emplace_back(Self());
    }

    template <typename F>
    void Map(F function) const {
        function(Get());
        for(auto& daughter : daughters) {
            daughter->Map(function);
        }
    }

    template <typename F>
    void Map_level(F function, size_t level=0) const {
        function(Get(), level);
        for(auto& daughter : daughters) {
            daughter->Map_level(function, level+1);
        }
    }

    /**
     * @brief Map_nodes runs through the tree recursively depth-first
     * @param function applied to each node
     * @note depth first is important for @c GetUniquePermutations
     */
    template <typename F>
    void Map_nodes(F function) const {
        for(auto& daughter : daughters) {
            daughter->Map_nodes(function);
        }
        function(Self());
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
                      [comp] (node_t d) { d->Sort(comp); });

        // then sort daughters itself
        daughters.sort([comp] (node_t a, node_t b) {
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
    bool IsEqual(const snode_t<U>& other, Compare comp) const {
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

    template<typename U>
    bool IsEqual(const snode_t<U>& other) const {
        return IsEqual(other, [] (const T& a, const U& b) { return a==b; });
    }

    /**
     * @brief DeepCopy uses the given transfrom from node_t to U to deeply copy the tree
     * @param transform function to convert node to new type U
     */
    template<typename U = T, typename Transform = std::function<U(node_t)> >
    snode_t<U> DeepCopy(Transform transform = [] (const node_t& n) { return n->Get(); } ) const {
        auto r = Tree<U>::MakeNode(transform(Self()));
        for(const auto& daughter : daughters) {
            r->AddDaughter(daughter->DeepCopy<U>(transform));
        }
        return r;
    }

    /**
     * @brief GetUniquePermutations calculates the unique permutations of the (assumed identical) leaves
     * @param leaves the leaves of the tree for convenience, defines also the order of indices in perms
     * @param perms permutations as indices corresponding to vector leaves
     *
     * @note We assume that the tree is already sorted with std::less comparison,
     * anything else will lead to wrong results.
     *
     * @todo extend function such that non-identical leaves are supported, i.e. let user provide comparison
     */
    void GetUniquePermutations(std::vector<node_t>& leaves,
                               std::vector<std::vector<size_t>>& perms) const
    {

        // ensure we're sorted (let's hope by standard std::less)
        // and everything is cleared
        Map_nodes([] (const node_t& n) {
            if(!n->is_sorted)
                throw std::runtime_error("Tree is not sorted");
        });
        leaves.resize(0);
        perms.resize(0);

        // append to each node of this tree some bitfield
        struct wrapped_t {
            wrapped_t(const node_t& n) : Node(n) {}
            node_t Node;
            std::uint64_t Bitfield = 0;
            bool operator<(const wrapped_t& rhs) const {
                // std::less here...
                const bool a_less_b = Node->Get() < rhs.Node->Get();
                const bool b_less_a = rhs.Node->Get() < Node->Get();
                if(!a_less_b && !b_less_a) {
                    return Bitfield < rhs.Bitfield;
                }
                return a_less_b;
            }
            bool operator==(const wrapped_t& rhs) const {
                // ignore leaves in do { ... } loop
                // Bitfield is only used there, so use this as a flag
                if(Bitfield>0) {
                    // in do {...} loop, we also assume the Node does not change
                    //assert(Node->Get() == rhs.Node->Get());
                    //assert(Node->Daughters().size() == rhs.Node->Daughters().size());
                    if(Node->IsLeaf())
                        return true;
                }
                return !(*this < rhs) && !(rhs < *this);
            }
        };

        // prepare the wrapped version of this tree
        // and also the leaves and the current permutation array
        using wrapped_node_t = typename Tree<wrapped_t>::node_t;
        wrapped_node_t wrapped = DeepCopy<wrapped_t>([] (const node_t& n) { return n; });

        using perm_t = std::vector<size_t>;
        perm_t current_perm;
        wrapped_node_t prev_leave = nullptr;

        std::vector<wrapped_node_t> wrapped_leaves;
        wrapped->Map_nodes([&wrapped_leaves] (wrapped_node_t node) {
            if(node->IsLeaf())
                wrapped_leaves.push_back(node);
        });

        for(unsigned i=0;i<wrapped_leaves.size();i++) {
            auto& leave = wrapped_leaves[i];
            if(i > sizeof(std::uint64_t)*8-1)
                throw std::runtime_error("Too many leaves");
            // important here to keep the order of leaves
            leaves.push_back(leave->Get().Node);
            current_perm.push_back(i);

            // do not use == here
            if(prev_leave && !(prev_leave->Get() == leave->Get())) {
                throw std::runtime_error("Not all leaves are equal");
            }
            else
                prev_leave = leave;
        }

        // remember the already found trees corresponding to permutations
        std::list<wrapped_node_t> unique_trees;

        // loop over all permuations
        do {
            // apply the current permutation to the leaves
            for(unsigned i=0;i<wrapped_leaves.size();i++) {
                auto& perm_leave = wrapped_leaves[current_perm[i]];
                perm_leave->Get().Bitfield = 1 << i;
            }

            // calculate the sum of upper bitfields,
            // we rely on deep-first recursion of Map_nodes
            wrapped->Map_nodes([] (wrapped_node_t node) {
                if(node->IsLeaf())
                    return;
                node->Get().Bitfield = 0;
                for(const auto& d : node->Daughters())
                    node->Get().Bitfield += d->Get().Bitfield;
            });

            // check that summation has worked
            assert(wrapped->Get().Bitfield == (unsigned)(1 << wrapped_leaves.size())-1 );

            // sort by Bitfield sums (but keep ordering of original nodes)
            wrapped->Sort();

            // have a look if we know this tree already
            bool found = false;
            for(const wrapped_node_t& t : unique_trees) {
                if(t->IsEqual(wrapped)) {
                    found = true;
                    break;
                }
            }
            if(!found) {
                auto t = wrapped->DeepCopy();
                t->Sort();
                unique_trees.emplace_back(move(t));
                // add it to the returned perms as well
                perms.emplace_back(current_perm);
            }

        }
        while(std::next_permutation(current_perm.begin(), current_perm.end()));

    }

};

} // namespace ant
