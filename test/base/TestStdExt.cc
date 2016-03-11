#include "catch.hpp"
#include "catch_config.h"

#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "base/std_ext/vector.h"
#include "base/std_ext/system.h"
#include "base/std_ext/shared_ptr_container.h"

#include "base/tmpfile_t.h"

#include <string>
#include <memory>
#include <iostream>

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
void TestLsFiles();
void TestSharedPtrContainer();

TEST_CASE("make_unique", "[base/std_ext]") {
    TestMakeUnique();
}

TEST_CASE("string stuff", "[base/std_ext") {
    TestString();
}

TEST_CASE("vector stuff", "[base/std_ext]") {
    TestVector();
}

TEST_CASE("system lsFiles", "[base/std_ext]") {
    TestLsFiles();
}

TEST_CASE("shared_ptr_container", "[base/std_ext]") {
    TestSharedPtrContainer();
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



void TestLsFiles() {
    tmpfolder_t folder;
    tmpfile_t a(folder, ".tst");
    a.write_testdata();

    list<string> files;

    REQUIRE_NOTHROW(files = std_ext::system::lsFiles(folder.foldername,".tst"));
    REQUIRE(files.size()==1);
    REQUIRE(files.front()==string(a.filename));


    tmpfile_t b(folder, ".tst");
    b.write_testdata();

    REQUIRE_NOTHROW(files = std_ext::system::lsFiles(folder.foldername,".tst"));
    REQUIRE(files.size()==2);
}

struct int_t {
    static size_t n_constructed;
    int val;
    int_t(int v) : val(v) { n_constructed++; }
    operator int() const { return val; }
};

size_t int_t::n_constructed = 0;

void TestSharedPtrContainer() {

    std_ext::shared_ptr_container<int_t, std::list> c1;
    c1.emplace_back(5);
    REQUIRE(*c1.begin() == 5);
    REQUIRE(c1[0] == 5);
    c1.emplace_back(6);
    c1.erase(c1.begin());
    REQUIRE(c1[0]==6);

    std_ext::shared_ptr_container<int_t, std::vector> c2{c1.begin(), c1.begin()};
    REQUIRE(c2[0]==6);
    REQUIRE(c2[1]==6);

    // number of emplace_back calls!
    REQUIRE(int_t::n_constructed == 2);

}
