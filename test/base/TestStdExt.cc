#include "catch.hpp"
#include "catch_config.h"

#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "base/std_ext/vector.h"

#include <string>
#include <memory>

using namespace std;
using namespace ant;

struct MemtestDummy {
    static unsigned n;
    MemtestDummy() { ++n; }
    ~MemtestDummy() { --n; }
};

unsigned MemtestDummy::n = 0;

void TestMakeUnique();
void TestString();
void TestVector();

TEST_CASE("make_unique", "[base/stdext]") {
    TestMakeUnique();
}

TEST_CASE("string stuff", "[base/stdext") {
    TestString();
}

TEST_CASE("vector stuff", "[base/stdext]") {
    TestVector();
}


void TestMakeUnique() {
    std::unique_ptr<MemtestDummy> d;

    REQUIRE_NOTHROW(d = std_ext::make_unique<MemtestDummy>());
    REQUIRE(MemtestDummy::n == 1);
    REQUIRE_NOTHROW(d = nullptr);
    REQUIRE(MemtestDummy::n == 0);
}

void TestString() {

    // ends_with
    const string test("Hallo Welt!");
    REQUIRE(std_ext::string_ends_with(test, "Welt!"));
    REQUIRE_FALSE(std_ext::string_ends_with(test, "Bla"));

    // removesubstr
    string str("foo bar foo bar!");
    REQUIRE_NOTHROW(std_ext::removesubstr(str,"bar"));
    REQUIRE(str == string("foo  foo !"));

    // sanitize
    const char* str2  = "\t\tHallo\n";
    string s;
    REQUIRE_NOTHROW(s=std_ext::string_sanitize(str2));
    REQUIRE(s == string("Hallo"));

    // formatter
    REQUIRE_NOTHROW(s = std_ext::formatter() << "hallo" << 2 << 5 << "du " << setw(3) << 1);
    REQUIRE(s == string("hallo25du   1"));

}


void TestVector() {
    std::vector<int> t = {1,2,3,4,5,6,7};

    REQUIRE(std_ext::contains(t,3));
    REQUIRE_FALSE(std_ext::contains(t,8));
}
