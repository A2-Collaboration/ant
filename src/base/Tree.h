#pragma once

// ignore warnings from library
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include "cereal/access.hpp"
#pragma GCC diagnostic pop

#include "std_ext/variadic.h"

#include <memory>
#include <list>
#include <vector>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <type_traits>

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

    template<typename Daughter, typename... Daughters>
    static typename
    std::enable_if<std_ext::first_element_is_tuple<Daughters...>::value, void>::type
    make_impl(node_t& head, const std::tuple<Daughter, Daughters...>& daughters) {
        auto& d = make_impl(head, std::get<0>(daughters));
        auto t = std_ext::strip_first_from_tuple(daughters);
        make_impl(d, std::get<0>(t));
        make_impl(head, std_ext::strip_first_from_tuple(t));
    }

    template<typename Daughter, typename... Daughters>
    static typename
    std::enable_if<!std_ext::first_element_is_tuple<Daughters...>::value, void>::type
    make_impl(node_t& head, const std::tuple<Daughter, Daughters...>& daughters) {
        make_impl(head, std::get<0>(daughters));
        make_impl(head, std_ext::strip_first_from_tuple(daughters));
    }

    template<typename Daughter = void>
    static void
    make_impl(node_t&, const std::tuple<>&) {}

    template<typename Item>
    static typename
    std::enable_if<std::is_convertible<Item, T>::value, node_t&>::type
    make_impl(node_t& head, Item&& item) {
        return head->CreateDaughter(std::forward<Item>(item));
    }

    template<typename Item>
    static typename
    std::enable_if<std::is_convertible<Item, node_t>::value, node_t&>::type
    make_impl(node_t& head, Item&& item) {
        return head->AddDaughter(std::forward<Item>(item)->DeepCopy());
    }

public:


    template <typename ... args_t>
    static node_t MakeNode(args_t&&... args) {
        // expose protected ctor by inheritance
        struct make_shared_enabler : Tree {
            make_shared_enabler(T&& data_) :
                Tree(std::forward<T>(data_)) {}
        };
        auto n = std::make_shared<make_shared_enabler>(T(std::forward<args_t>(args)...));
        n->self = n;
        return n;
    }

    template<typename Head, typename... Daughters>
    static node_t
    Make(Head&& head, const std::tuple<Daughters...>& daughters) {
        auto headnode = MakeNode(std::forward<Head>(head));
        make_impl(headnode, daughters);
        return headnode;
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

    /**
     * @brief Unlink remove this node from parent
     * @note ensure that the instance of the node is local, and NOT some
     * referenced node in the original tree
     */
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
    template<typename U = T, typename Transform = std::function<U(const node_t&)> >
    snode_t<U> DeepCopy(Transform transform) const {
        auto r = Tree<U>::MakeNode(transform(Self()));
        for(const auto& daughter : daughters) {
            r->AddDaughter(daughter->template DeepCopy<U>(transform));
        }
        return r;
    }

    template<typename U = T>
    snode_t<U> DeepCopy() const {
        auto r = Tree<U>::MakeNode(Get());
        for(const auto& daughter : daughters) {
            r->AddDaughter(daughter->template DeepCopy<U>());
        }
        return r;
    }

    /**
     * @brief GetUniquePermutations calculates the unique permutations of the leaves
     * @param leaves the leaves of the tree for convenience, defines also the order of indices in perms
     * @param perms permutations as indices corresponding to vector leaves (offseted by possible i_leave_offset).
     * Use indices to assign particles to leaves (not leaves to particles)
     * i_leave_offset separates the trivial equivalence classes from the actually permuted leaves
     *
     * @note We assume that the tree is already sorted with std::less comparison,
     * anything else will lead to wrong results.
     *
     */
    void GetUniquePermutations(std::vector<node_t>& leaves,
                               std::vector<std::vector<int>>& perms,
                               int& i_leave_offset) const
    {

        // ensure we're sorted (let's hope by standard std::less of underyling type)
        Map_nodes([] (const node_t& n) {
            if(!n->is_sorted)
                throw std::runtime_error("Tree is not sorted");
        });

        // be nice and convenient :)
        leaves.resize(0);
        perms.resize(0);

        // append to each node of this tree some bitfield
        struct wrapped_t {
            wrapped_t(const node_t& n) : Node(n) {}
            node_t Node;
            std::uint64_t Bitfield = 0;
        };

        // prepare the wrapped version of this tree
        // and also the leaves and the current permutation array
        using wrapped_node_t = typename Tree<wrapped_t>::node_t;
        wrapped_node_t wrapped = DeepCopy<wrapped_t>([] (const node_t& n) { return n; });

        // figure out how many equivalence classes we have in the leaves
        // and how many representatives are there
        // the following algorithm in the do {...} loop
        // only supports one equivalence class with more than one member

        const auto equi_relation = [] (const wrapped_node_t& a, const wrapped_node_t& b) {
            const auto& a_ = a->Get().Node->Get();
            const auto& b_ = b->Get().Node->Get();
            return !(a_ < b_) && !(b_ < a_);
        };

        using equi_class_t = std::vector<wrapped_node_t>;
        std::vector<equi_class_t> equi_classes;

        unsigned n_leaves = 0;

        wrapped->Map_nodes(
                    [&equi_classes, &n_leaves, equi_relation] (const wrapped_node_t& node) {
            // only leaves
            if(!node->IsLeaf())
                return;
            n_leaves++;
            for(auto& equi : equi_classes) {
                // compare with underlying node
                if(equi_relation(equi.front(), node)) {
                    equi.emplace_back(node);
                    return;
                }
            }
            // not found in equivalence_classes, add new class with one entry
            equi_classes.emplace_back(equi_class_t{node});
        });

        assert(!equi_classes.empty());

        if(n_leaves > sizeof(std::uint64_t)*8)
            throw std::runtime_error("More leaves than 64bit field can cope with");

        // sort equivalence classes by number of members, highest last
        std::sort(equi_classes.begin(), equi_classes.end(),
                  [] (const equi_class_t& a, const equi_class_t& b) {
            // make sure sort is stable
            if(a.size() == b.size()) {
                return a.front()->Get().Node->Get() < b.front()->Get().Node->Get();
            }
            return a.size() < b.size();
        });

        // check number of equi classes with size() == 1
        auto n_size_1 = std::count_if(equi_classes.begin(), equi_classes.end(),
                                      [] (const equi_class_t& e) {
            return e.size() == 1;
        });

        if(equi_classes.size() > 1 && n_size_1 != int(equi_classes.size())-1)
            throw std::runtime_error("Found leave equivalence classes with too many different types");

        // fill the leaves, important to keep order as defined by equi_classes
        // other the offset does not make sense!
        i_leave_offset = equi_classes.size()-1;
        for(const auto& equi : equi_classes) {
            for(const auto& leave : equi)
                leaves.emplace_back(leave->Get().Node);
        }

        // from now on, only work on the non-trivial equiclass,
        // called wrapped_leaves
        // and create the trivial permutation
        const auto& wrapped_leaves = equi_classes.back();
        using perm_t = std::vector<int>;
        perm_t current_perm(wrapped_leaves.size());
        std::iota(current_perm.begin(), current_perm.end(), 0);

        // the following loop over the permutations
        // needs those sorting and equality relations
        // depending on the Bitfield
        const auto wrapped_less = [] (const wrapped_t& a, const wrapped_t& b) {
            // std::less here...
            const bool a_less_b = a.Node->Get() < b.Node->Get();
            const bool b_less_a = b.Node->Get() < a.Node->Get();
            if(!a_less_b && !b_less_a) {
                return a.Bitfield < b.Bitfield;
            }
            return a_less_b;
        };

        const auto wrapped_equal = [wrapped_less] (const wrapped_t& a, const wrapped_t& b) {
            //assert(Node->Get() == rhs.Node->Get());
            //assert(Node->Daughters().size() == rhs.Node->Daughters().size());
            // ignore leaves in do { ... } loop
            if(a.Node->IsLeaf() || b.Node->IsLeaf())
                return true;
            return !wrapped_less(a, b) && !wrapped_less(b, a);
        };

        // remember the already found trees corresponding to permutations
        // uniqueness is defined by wrapped_equal and eventually
        // by bitfield sums
        std::list<wrapped_node_t> unique_trees;

        // loop over all permuations
        do {
            // apply the current permutation to the leaves
            for(unsigned i=0;i<wrapped_leaves.size();i++) {
                auto& leave = wrapped_leaves[i]->Get();
                leave.Bitfield = 1 << current_perm[i];
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
            wrapped->Sort(wrapped_less);

            // have a look if we know this tree already
            bool found = false;
            for(const wrapped_node_t& t : unique_trees) {
                if(t->IsEqual(wrapped, wrapped_equal)) {
                    found = true;
                    break;
                }
            }
            if(!found) {
                auto t = wrapped->DeepCopy();
                t->Sort(wrapped_less);
                unique_trees.emplace_back(move(t));
                // add it to the returned perms as well
                perms.emplace_back(current_perm);
            }

        }
        while(std::next_permutation(current_perm.begin(), current_perm.end()));
    }

};

} // namespace ant
