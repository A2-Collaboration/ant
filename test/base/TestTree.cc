#include "catch.hpp"
#include "base/Tree.h"
#include <iostream>

using namespace std;
using namespace ant;

TEST_CASE("Tree: Default ctor", "[base]") {
//    REQUIRE_NOTHROW( auto a = Tree<int>::makeNode(); );
    REQUIRE_NOTHROW( auto a = Tree<int>::makeNode(10); );
}

TEST_CASE("Tree: Assignment", "[base]") {
    auto a = Tree<int>::makeNode(10);

    REQUIRE(**a == 10);
    REQUIRE(a->get() == 10);
    REQUIRE_NOTHROW(**a = 20;);
    REQUIRE(**a == 20);
}

TEST_CASE("Tree: Parents", "[base]") {
    auto a = Tree<int>::makeNode(10);
    REQUIRE(a->isLeaf());
    REQUIRE(a->isRoot());

    REQUIRE_NOTHROW(auto b = a->createDaughter(30);
    REQUIRE_FALSE(a->isLeaf());
    REQUIRE(a->isRoot());

    REQUIRE_FALSE(b->isRoot());
    REQUIRE(b->isLeaf());

    auto c = a->createDaughter(30);
    REQUIRE(**c==30);
    REQUIRE_FALSE(c->isRoot());
    REQUIRE(c->isLeaf());
    REQUIRE(a->Daughters().size()==2);
    );

}

TEST_CASE("Tree: Unlink", "[base]") {
    auto a = Tree<int>::makeNode(10);
    auto b = a->createDaughter(20);
    auto c = b->createDaughter(30);
    auto d = b->createDaughter(35);
    REQUIRE(a->isRoot());
    REQUIRE_FALSE(a->isLeaf());
    REQUIRE_FALSE(b->isLeaf());
    REQUIRE_FALSE(b->isRoot());
    b->Unlink();
    REQUIRE(b->isRoot());
    REQUIRE(a->isLeaf());
}

TEST_CASE("Tree: map", "[base]") {
    auto a = Tree<int>::makeNode(10);
    auto b = a->createDaughter(20);
    auto c = b->createDaughter(30);
    auto d = b->createDaughter(35);
    int sum = 0;
    a->map([&sum] (int i) { sum += i;});
    REQUIRE(sum == 95);
}

TEST_CASE("Tree: size", "[base]") {
    auto a = Tree<int>::makeNode(10);
    auto b = a->createDaughter(20);
    auto c = b->createDaughter(30);
    auto d = b->createDaughter(35);
    REQUIRE(a->size() == 4);
}

TEST_CASE("Tree: depth", "[base]") {
    auto a = Tree<int>::makeNode(10);
    auto b = a->createDaughter(20);
    auto c = b->createDaughter(30);
    auto d = b->createDaughter(35);
    REQUIRE(a->depth() == 3);
}

TEST_CASE("Tree: sort", "[base]") {
    auto int_sorter = [] (int a, int b) { return a<b; };

    auto a = Tree<int>::makeNode(1);
    auto a0 = a->createDaughter(2);
    auto a1 = a->createDaughter(1);

    auto a00 = a0->createDaughter(2);
    auto a01 = a0->createDaughter(2);

    auto a10 = a1->createDaughter(3);
    auto a11 = a1->createDaughter(3);

    a->sort(int_sorter);

    auto a0_ = a->Daughters().front();
    REQUIRE(**a0_ == 1);

    auto a01_ = a->Daughters().front()->Daughters().back();
    REQUIRE(**a01_ == 3);


    auto b = Tree<int>::makeNode(1);
    auto b0 = b->createDaughter(1);
    auto b1 = b->createDaughter(1);

    auto b00 = b0->createDaughter(3);
    auto b01 = b0->createDaughter(3);

    auto b10 = b1->createDaughter(2);
    auto b11 = b1->createDaughter(2);

    b->sort(int_sorter);

    auto b0_ = b->Daughters().front();
    REQUIRE(**b0_ == 1);

    auto b00_ = b->Daughters().front()->Daughters().front();
    REQUIRE(**b00_ == 2);
}

