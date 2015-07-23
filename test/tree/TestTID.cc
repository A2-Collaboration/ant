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
    TID a(0,0);
    TID b(0,1);
    TID c(0,0,true);

    REQUIRE(a == a);

    REQUIRE_FALSE(a == b);
    REQUIRE_FALSE( a == c);

    REQUIRE(a != b);
    REQUIRE(a != c);

    REQUIRE(a < b);
    REQUIRE(b > a);
    REQUIRE(a < c);
    REQUIRE(b < c);
}
