#include "catch.hpp"
#include "base/BinSettings.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
#include <iostream>
#include <sstream>

using namespace std;
using namespace ant;

TEST_CASE("BinSettings: Default ctor", "[base]") {
    REQUIRE_NOTHROW( std_ext::make_unique<BinSettings>(20,0,10) );
    REQUIRE_NOTHROW( std_ext::make_unique<BinSettings>(20) );
}

TEST_CASE("BinSettings: Parse from string", "[base]") {

    const BinSettings a(10, -5, 8.4);
    const BinSettings c(2, 1, 3);
    stringstream ss_a;
    REQUIRE_NOTHROW(ss_a << a);

    ss_a.seekg(0);
    BinSettings b(0);
    REQUIRE(ss_a >> b);
    REQUIRE(b == a);
    REQUIRE(c != a);

   stringstream ss_b;
   ss_b << "(200,[0.95,1.2])";ss_b.seekg(0);
   REQUIRE(ss_b >> b);
   REQUIRE(b.Bins()==200);
}
