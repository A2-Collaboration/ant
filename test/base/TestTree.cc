#include "catch.hpp"
#include "base/Tree.h"
#include <iostream>

#include "base/ParticleTypeTree.h"
#include "analysis/utils/ParticleTools.h"

using namespace std;
using namespace ant;

TEST_CASE("Tree: Default ctor", "[base]") {
    REQUIRE_NOTHROW(Tree<int>::MakeNode(10) );

    struct Object_t {
        int A;
        int B;
        Object_t(int a, int b) : A(a), B(b) {}
    };
    REQUIRE_NOTHROW(Tree<Object_t>::MakeNode(3,4));
}

TEST_CASE("Tree: Make from nested tuples", "[base]") {
    {
        auto t = Tree<int>::Make(1, std::make_tuple(2, 4, 3));
        CHECK(t->Daughters().size() == 3);
        CHECK(t->Get() == 1);
    }
    {
        auto t = Tree<int>::Make(1, std::make_tuple(4, std::make_tuple(2, 17), 3));
        CHECK(t->Daughters().size() == 2);
        CHECK(t->Daughters().front()->Daughters().size() == 2);
        CHECK(t->Daughters().front()->Daughters().back()->Get() == 17);
    }
    {
        using PTree_t = ParticleTypeTree::element_type;

        auto& gp  = ParticleTypeDatabase::BeamProton;
        auto  g   = PTree_t::MakeNode(ParticleTypeDatabase::Photon);
        auto  p   = PTree_t::MakeNode(ParticleTypeDatabase::Proton);
        auto  pi0 = PTree_t::MakeNode(ParticleTypeDatabase::Pi0);

        auto t = PTree_t::Make(
                     gp,
                     std::make_tuple(p, pi0, std::make_tuple(g, g))
                     );
        t->Sort();
        // currently assumes the trees in database are build differently
        CHECK(t->IsEqual(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0_2g)));
    }
}




TEST_CASE("Tree: Assignment", "[base]") {
    auto a = Tree<int>::MakeNode(10);

    REQUIRE(**a == 10);
    REQUIRE(a->Get() == 10);
    **a = 20;
    REQUIRE(**a == 20);
}

TEST_CASE("Tree: Parents", "[base]") {
    auto a = Tree<int>::MakeNode(10);
    REQUIRE(a->IsLeaf());
    REQUIRE(a->IsRoot());

    auto b = a->CreateDaughter(30);
    REQUIRE_FALSE(a->IsLeaf());
    REQUIRE(a->IsRoot());

    REQUIRE_FALSE(b->IsRoot());
    REQUIRE(b->IsLeaf());

    auto c = a->CreateDaughter(30);
    REQUIRE(**c==30);
    REQUIRE_FALSE(c->IsRoot());
    REQUIRE(c->IsLeaf());
    REQUIRE(a->Daughters().size()==2);

}

TEST_CASE("Tree: Unlink", "[base]") {
    auto a = Tree<int>::MakeNode(10);
    auto b = a->CreateDaughter(20);
    auto c = b->CreateDaughter(30);
    auto d = b->CreateDaughter(35);
    REQUIRE(a->IsRoot());
    REQUIRE_FALSE(a->IsLeaf());
    REQUIRE_FALSE(b->IsLeaf());
    REQUIRE_FALSE(b->IsRoot());
    b->Unlink();
    REQUIRE(b->IsRoot());
    REQUIRE(a->IsLeaf());
}

TEST_CASE("Tree: map", "[base]") {
    auto a = Tree<int>::MakeNode(10);
    auto b = a->CreateDaughter(20);
    auto c = b->CreateDaughter(30);
    auto d = b->CreateDaughter(35);
    int sum = 0;
    a->Map([&sum] (int i) { sum += i;});
    REQUIRE(sum == 95);
}

TEST_CASE("Tree: size", "[base]") {
    auto a = Tree<int>::MakeNode(10);
    auto b = a->CreateDaughter(20);
    auto c = b->CreateDaughter(30);
    auto d = b->CreateDaughter(35);
    REQUIRE(a->Size() == 4);
}

TEST_CASE("Tree: depth", "[base]") {
    auto a = Tree<int>::MakeNode(10);
    auto b = a->CreateDaughter(20);
    auto c = b->CreateDaughter(30);
    auto d = b->CreateDaughter(35);
    REQUIRE(a->Depth() == 3);
}

TEST_CASE("Tree: sort", "[base]") {
    auto int_sorter = [] (int a, int b) { return a<b; };

    auto a = Tree<int>::MakeNode(1);
    auto a0 = a->CreateDaughter(2);
    auto a1 = a->CreateDaughter(1);

    auto a00 = a0->CreateDaughter(2);
    auto a01 = a0->CreateDaughter(2);

    auto a10 = a1->CreateDaughter(3);
    auto a11 = a1->CreateDaughter(3);

    a->Sort(int_sorter);

    auto a0_ = a->Daughters().front();
    REQUIRE(**a0_ == 1);

    auto a01_ = a->Daughters().front()->Daughters().back();
    REQUIRE(**a01_ == 3);


    auto b = Tree<int>::MakeNode(1);
    auto b0 = b->CreateDaughter(1);
    auto b1 = b->CreateDaughter(1);

    auto b00 = b0->CreateDaughter(3);
    auto b01 = b0->CreateDaughter(3);

    auto b10 = b1->CreateDaughter(2);
    auto b11 = b1->CreateDaughter(2);

    b->Sort(int_sorter);

    auto b0_ = b->Daughters().front();
    REQUIRE(**b0_ == 1);

    auto b00_ = b->Daughters().front()->Daughters().front();
    REQUIRE(**b00_ == 2);
}

TEST_CASE("Tree: compare", "[base]") {
    auto a = Tree<int>::MakeNode(1);
    auto a0 = a->CreateDaughter(1);
    auto a1 = a->CreateDaughter(1);
    a0->CreateDaughter(2);
    a0->CreateDaughter(2);
    a1->CreateDaughter(3);
    a1->CreateDaughter(3);

    auto b = Tree<int>::MakeNode(1);
    auto b0 = b->CreateDaughter(1);
    auto b1 = b->CreateDaughter(1);
    b0->CreateDaughter(3);
    b0->CreateDaughter(3);
    b1->CreateDaughter(2);
    b1->CreateDaughter(2);

    auto int_comparer = [] (int a, int b) { return a==b; };
    auto int_sorter = [] (int a, int b) { return a<b; };


    REQUIRE_THROWS(a->IsEqual(b, int_comparer));
    a->Sort(int_sorter);
    REQUIRE_THROWS(a->IsEqual(b, int_comparer));
    b->Sort(int_sorter);
    REQUIRE_NOTHROW(a->IsEqual(b, int_comparer));


    REQUIRE(a->IsEqual(b, int_comparer));

    // check the number of calls of the comparer
    size_t count = 0;
    auto int_comparer_count = [&count] (int a, int b) {
        count++;
        return a == b;
    };
    REQUIRE(a->IsEqual(b, int_comparer_count));
    REQUIRE(count == a->Size());


    // modify the trees, check if sorted flag is unset
    // which makes Compare() throw an exception

    a0->CreateDaughter(1);
    b1->CreateDaughter(1);
    REQUIRE_THROWS(a->IsEqual(b, int_comparer));
    a->Sort(int_sorter);
    b->Sort(int_sorter);
    REQUIRE_NOTHROW(a->IsEqual(b, int_comparer));
    REQUIRE(a->IsEqual(b, int_comparer));


    // make the trees unequal
    a0->CreateDaughter(1);
    b1->CreateDaughter(2);
    a->Sort(int_sorter);
    b->Sort(int_sorter);
    REQUIRE_FALSE(a->IsEqual(b, int_comparer));
}

TEST_CASE("Tree: DeepCopy", "[base]") {
    auto a = Tree<int>::MakeNode(1);
    auto a0 = a->CreateDaughter(1);
    auto a1 = a->CreateDaughter(1);
    a0->CreateDaughter(2);
    a0->CreateDaughter(2);
    a1->CreateDaughter(3);
    a1->CreateDaughter(3);

    auto copy = a->DeepCopy<double>();

    REQUIRE(copy->Get() == 1);
    REQUIRE(copy->Daughters().front()->Daughters().front()->Get() == 2);
    REQUIRE(copy->Daughters().back()->Daughters().front()->Get() == 3);

}

void dotest_perms();

TEST_CASE("Tree: GetUniquePermutations", "[base]") {
    dotest_perms();
}

void dotest_perms() {
    using perms_t = std::vector<std::vector<int>>;
    perms_t perms;

    auto check_perms = [&perms] (const perms_t& expected) {
        auto it1 = perms.begin();
        auto it2 = expected.begin();
        while(it1 != perms.end() && it2 != expected.end()) {
            REQUIRE(*it1 == *it2);
            ++it1;
            ++it2;
        }
    };

    vector<ParticleTypeTree> leaves_ptree;
    ParticleTypeTree ptree;
    int i_leave_offset;

    ptree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_2Pi0Eta_6g);
    REQUIRE_NOTHROW(ptree->GetUniquePermutations(leaves_ptree, perms, i_leave_offset));
    REQUIRE(i_leave_offset == 1);
    REQUIRE(leaves_ptree.size() == 7);
    REQUIRE(perms.size() == 45);

    ptree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_3Pi0_6g);
    REQUIRE_NOTHROW(ptree->GetUniquePermutations(leaves_ptree, perms, i_leave_offset));
    REQUIRE(i_leave_offset == 1);
    REQUIRE(leaves_ptree.size() == 7);
    REQUIRE(perms.size() == 15);

    ptree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::EtaPrime_gOmega_ggPi0_4g);
    REQUIRE_NOTHROW(ptree->GetUniquePermutations(leaves_ptree, perms, i_leave_offset));
    REQUIRE(i_leave_offset == 1);
    REQUIRE(leaves_ptree.size() == 5);
    REQUIRE(perms.size() == 12);

    ptree = ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::SigmaPlusK0s_6g);
    REQUIRE_NOTHROW(ptree->GetUniquePermutations(leaves_ptree, perms, i_leave_offset));
    REQUIRE(i_leave_offset == 1);
    REQUIRE(leaves_ptree.size() == 7);
    REQUIRE(perms.size() == 45);


    auto a = Tree<int>::MakeNode(0);
    auto a0 = a->CreateDaughter(1);
    auto a1 = a->CreateDaughter(1);
    a0->CreateDaughter(2);
    a0->CreateDaughter(2);
    a1->CreateDaughter(2);
    a1->CreateDaughter(2);
    a->CreateDaughter(3);
    a->CreateDaughter(4);
    a->CreateDaughter(5);

    vector<Tree<int>::node_t> leaves_int;

    a->Sort();
    REQUIRE_NOTHROW(a->GetUniquePermutations(leaves_int, perms, i_leave_offset));

    REQUIRE(leaves_int.size() == 7);
    REQUIRE(perms.size() == 3);
    REQUIRE(i_leave_offset == 3);
    check_perms({ {0,1,2,3}, {0,2,1,3}, {0,3,1,2} });


    // test some unsupported tree (too many leaves with different types aka numbers)
    auto b = Tree<int>::MakeNode(0);
    b->CreateDaughter(1);
    b->CreateDaughter(1);
    b->CreateDaughter(2);
    b->CreateDaughter(3);
    b->CreateDaughter(3);

    REQUIRE_THROWS(b->GetUniquePermutations(leaves_int, perms, i_leave_offset));


    // simple struct with just std::less as operator (no operator==)
    struct my_t {
        int a;
        my_t(int a_) : a(a_) {}
        bool operator<(const my_t& b) const {
            return a < b.a;
        }
    };

    auto c = Tree<my_t>::MakeNode(0);
    c->CreateDaughter(1);
    c->CreateDaughter(2);
    c->CreateDaughter(2);
    c->Sort();

    vector<Tree<my_t>::node_t> leaves_my_t;
    REQUIRE_NOTHROW(c->GetUniquePermutations(leaves_my_t, perms, i_leave_offset));
    REQUIRE(perms.size() == 1);
}
