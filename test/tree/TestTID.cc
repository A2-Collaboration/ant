#include "catch.hpp"

#include "TDataRecord.h"

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("TID","[operators]")
{
    dotest();
}

void dotest()
{
    TID a(0);
    TID b(1);
    TID c(0,true);

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
    REQUIRE(a.Value == 1);
    REQUIRE(b.Value == 0);

    TID d;

    // copy
    REQUIRE_NOTHROW(d = c);

    // flags same after copy
    REQUIRE_NOTHROW(d.Flags = c.Flags);

    //flags same after increment
    REQUIRE((++c).Flags == d.Flags);
}
