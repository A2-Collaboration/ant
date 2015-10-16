#include "catch.hpp"
#include "base/Tree.h"
#include <iostream>

using namespace std;
using namespace ant;

TEST_CASE("Tree: Default ctor", "[base]") {
    REQUIRE_NOTHROW( auto a = Tree<int>::MakeNode(10); );

    struct Object_t {
        int A;
        int B;
        Object_t(int a, int b) : A(a), B(b) {}
    };
    REQUIRE_NOTHROW(auto b = Tree<Object_t>::MakeNode(3,4));
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
