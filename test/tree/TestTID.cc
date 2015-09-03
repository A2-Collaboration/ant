#include "catch.hpp"

#include "TDataRecord.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("TID","[operators]")
{
    dotest();
}

void dotest()
{
    TID a(0, 0u);
    TID b(0, 1u);
    TID c(0, 0,true);

    cout << a << endl;
    cout << b << endl;
    cout << c << endl;

    REQUIRE(a == a);

    REQUIRE_FALSE(a == b);
    REQUIRE_FALSE(a == c);

    REQUIRE(a != b);
    REQUIRE(a != c);

    REQUIRE(a < b);
    REQUIRE(b > a);
    REQUIRE(a < c);
    REQUIRE(b < c);

    // increment / decrement
    REQUIRE(++a == b);
    REQUIRE(--b != a);
    REQUIRE(a.Lower == 1);
    REQUIRE(b.Lower == 0);

    TID d;

    // copy
    REQUIRE_NOTHROW(d = c);

    // flags same after copy
    REQUIRE_NOTHROW(d.Flags = c.Flags);

    //flags same after increment
    REQUIRE((++c).Flags == d.Flags);

    // invalid
    TID invalid;
    REQUIRE(invalid.IsInvalid());
    REQUIRE(invalid != invalid);
    REQUIRE_FALSE(invalid < invalid);
    REQUIRE_FALSE(invalid > invalid);
    REQUIRE_FALSE(invalid <= invalid);
    REQUIRE_FALSE(invalid >= invalid);
    REQUIRE_FALSE(a == invalid);
    REQUIRE_FALSE(invalid == a);
    REQUIRE_FALSE(a < invalid);
    REQUIRE_FALSE(a > invalid);
    REQUIRE_FALSE(invalid < a);
    REQUIRE_FALSE(invalid > a);
    REQUIRE_FALSE(a <= invalid);
    REQUIRE_FALSE(a >= invalid);
    REQUIRE_FALSE(invalid <= a);
    REQUIRE_FALSE(invalid >= a);
}
