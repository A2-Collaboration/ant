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

    auto b = Tree<int>::makeNode(20);
    REQUIRE_NOTHROW(Tree<int>::link(a,b));
    REQUIRE(a->isLeaf() == false);
    REQUIRE(a->isRoot() == true);

    REQUIRE(b->isRoot()==false);
    REQUIRE(b->isLeaf()==true);

    auto c = Tree<int>::makeNode(a,30);
    REQUIRE(**c==30);
    REQUIRE(c->isRoot() == false);
    REQUIRE(c->isLeaf() == true);
    REQUIRE(a->Daughters().size()==2);

}
