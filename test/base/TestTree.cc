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
    REQUIRE(a->isLeaf() == true);
    REQUIRE(a->isRoot() == true);

    REQUIRE_NOTHROW(auto b = a->createDaughter(30);
    REQUIRE(a->isLeaf() == false);
    REQUIRE(a->isRoot() == true);

    REQUIRE(b->isRoot()==false);
    REQUIRE(b->isLeaf()==true);

    auto c = a->createDaughter(30);
    REQUIRE(**c==30);
    REQUIRE(c->isRoot() == false);
    REQUIRE(c->isLeaf() == true);
    REQUIRE(a->Daughters().size()==2);
    );

}

TEST_CASE("Tree: Unlink", "[base]") {
    auto a = Tree<int>::makeNode(10);
    auto b = a->createDaughter(20);
    auto c = b->createDaughter(30);
    auto d = b->createDaughter(35);
    REQUIRE(a->isRoot() == true);
    REQUIRE(a->isLeaf() == false);
    REQUIRE(b->isLeaf() == false);
    REQUIRE(b->isRoot() == false);
    b->Unlink();
    REQUIRE(b->isRoot() == true);
    REQUIRE(a->isLeaf() == true);
}

TEST_CASE("Tree: map", "[base]") {
    auto a = Tree<int>::makeNode(10);
    auto b = a->createDaughter(20);
    auto c = b->createDaughter(30);
    auto d = b->createDaughter(35);
    a->map([] (const int& i) { cout << i << endl;});
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
