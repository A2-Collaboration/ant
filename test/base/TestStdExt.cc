#include "catch.hpp"
#include "catch_config.h"

#include "base/std_ext/memory.h"
#include "base/std_ext/string.h"
#include "base/std_ext/container.h"
#include "base/std_ext/system.h"
#include "base/std_ext/shared_ptr_container.h"
#include "base/std_ext/math.h"
#include "base/std_ext/misc.h"
#include "base/std_ext/vector.h"
#include "base/std_ext/map.h"

#include "base/tmpfile_t.h"

#include <string>
#include <memory>
#include <iostream>
#include <random>
#include <map>

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
void TestMap();
void TestLsFiles();
void TestSharedPtrContainer();
void TestRMSIQR();
void TestDereference();

TEST_CASE("make_unique", "[base/std_ext]") {
    TestMakeUnique();
}

TEST_CASE("string stuff", "[base/std_ext") {
    TestString();
}

TEST_CASE("vector stuff", "[base/std_ext]") {
    TestVector();
}

TEST_CASE("map stuff", "[base/std_ext]") {
    TestMap();
}

TEST_CASE("system lsFiles", "[base/std_ext]") {
    TestLsFiles();
}

TEST_CASE("shared_ptr_container", "[base/std_ext]") {
    TestSharedPtrContainer();
}

TEST_CASE("RMS_t and IQR_t", "[base/std_ext]") {
    TestRMSIQR();
}

TEST_CASE("Dereference", "[base/std_ext]") {
    TestDereference();
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
    // starts_with
    REQUIRE(std_ext::string_starts_with(test, "Hallo"));

    // tokenize
    string str("foo bar foo bar!");
    const auto tokens = std_ext::tokenize_string(str, " ");
    REQUIRE(tokens.size() == 4);
    REQUIRE(tokens.back() == "bar!");
    // concantenate
    REQUIRE(std_ext::concatenate_string(tokens, "_") == "foo_bar_foo_bar!");
    // removesubstr
    REQUIRE_NOTHROW(std_ext::removesubstr(str,"bar"));
    REQUIRE(str == string("foo  foo !"));

    // sanitize
    const char* str2  = "\t\tHallo\n";
    string s;
    REQUIRE_NOTHROW(s=std_ext::string_sanitize(str2));
    REQUIRE(s == string("Hallo"));

    REQUIRE(std_ext::string_sanitize("   ") == "");
    REQUIRE(std_ext::string_sanitize(" b  \n\t") == "b");
    REQUIRE(std_ext::string_sanitize("abc") == "abc");
    REQUIRE(std_ext::string_sanitize("") == "");
    REQUIRE(std_ext::string_sanitize("a") == "a");
    REQUIRE(std_ext::string_sanitize("a b") == "a b");

    // remove characters
    string special_chars("test #pi^{0}");
    std_ext::remove_chars(special_chars, {'#', '{', '}', '^'});
    REQUIRE(special_chars == string("test pi0"));

    // formatter
    REQUIRE_NOTHROW(s = std_ext::formatter() << "hallo" << 2 << 5 << "du " << setw(3) << 1);
    REQUIRE(s == string("hallo25du   1"));

    // basename
    REQUIRE(std_ext::basename("just/some/path") == "/path");
    REQUIRE(std_ext::basename("no_path") == "no_path");

}


void TestVector() {
    std::vector<int> t = {1,2,3,4,5,6,7};

    REQUIRE(std_ext::contains(t,3));
    REQUIRE_FALSE(std_ext::contains(t,8));

    const std::vector<int> v = {3,6,2,8,4};
    const auto sorted_asc = std_ext::get_sorted_indices(v);
    const auto sorted_desc = std_ext::get_sorted_indices_desc(v);
    REQUIRE(sorted_asc.size() == 5);
    REQUIRE(sorted_asc.at(1) == 0);
    REQUIRE(sorted_asc[2] == 4);
    REQUIRE(v[sorted_asc.front()] == 2);
    REQUIRE(v[sorted_asc.back()] == 8);
    REQUIRE(sorted_desc.front() == 3);
    REQUIRE(sorted_desc.back() == 2);

    // test circular rotations / shifts in container.h
    auto s(v);
    std_ext::shift_right(s);
    REQUIRE(s[0] == 4);
    std_ext::shift_right(s);
    REQUIRE(s[0] == 8);
    std::list<int> l;
    std::copy(s.begin(), s.end(), std::back_inserter(l));
    std_ext::shift_left(l);
    REQUIRE(l.back() == 8);
    std_ext::shift_left(l);
    REQUIRE(l.front() == 3);
}


void TestMap() {
    std::map<int, double> m = {{1, .7},
                               {2, .4},
                               {3, .9},
                               {4,-.2},
                               {5, .3}};

    REQUIRE(std_ext::min_map_element(m)->first == 4);
    REQUIRE(std_ext::max_map_element(m)->first == 3);

    std::vector<double> v;
    std::transform(m.begin(), m.end(), std::back_inserter(v), std_ext::second(m));

    REQUIRE(v.front() == .7);
    REQUIRE(v.at(2)   == .9);
    REQUIRE(v.back()  == .3);
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
    int_t& operator++() {
        ++val; return *this;
    }
    int get_value() const { return *this; }
    explicit operator bool() const { return val == 0; }
};

size_t int_t::n_constructed = 0;

void TestSharedPtrContainer() {

    std_ext::shared_ptr_container<int_t, std::list> c1;
    c1.emplace_back(5);
    REQUIRE(*c1.begin() == 5);
    REQUIRE(c1.front() == 5);
    c1.emplace_back(6);
    c1.erase(c1.begin());
    REQUIRE(c1.front()==6);
    c1.emplace_back(7);

    // c1 == {6, 7}
    REQUIRE(c1.size()  == 2);
    REQUIRE(c1.front() == 6);
    REQUIRE(c1.back()  == 7);

    std_ext::shared_ptr_container<int_t, std::vector> c2;
    c2.insert(c2.end(), c1.begin(), c1.end());
    c2.insert(c2.end(), c1.begin(), c1.end());
    REQUIRE(c2.size()==4);

    for(auto& item : c2)
        ++item;
    // due to pointers, the elements 0,1 are identical to 2,3
    REQUIRE(c2.get_ptr_at(0));
    REQUIRE(c2.get_ptr_at(0) == c2.get_ptr_at(2));
    REQUIRE(c2.get_ptr_at(1) == c2.get_ptr_at(3));
    REQUIRE(c2[0]==8);
    REQUIRE(c2[1]==9);

    auto ptr = c1.begin().get_ptr();
    using ptr_t = decltype(ptr);
    REQUIRE(std::is_const<ptr_t::element_type>::value);
    auto is_convertible = std::is_convertible<ptr_t,shared_ptr<int_t>>::value;
    REQUIRE(!is_convertible);

    // use range construct
    const std_ext::shared_ptr_container<int_t, std::list> c3(c2.begin(), c2.end());
    REQUIRE(c3.size()==4);

    // make sure the constness of the container is correctly propagated to the iterator
    constexpr auto is_assignable = std::is_assignable<decltype(*c3.begin()),int_t>::value;
    REQUIRE_FALSE(is_assignable);

    unsigned n = 0;
    for(auto i : c3.get_iter()) {
        REQUIRE(i.get_ptr());
        ++n;
    }
    REQUIRE(n == 4);

    unsigned n_eights = 0;
    for(auto i : c3.get_iter([] (int_t i) { return i==8;} )) {
        REQUIRE(i.get_ptr());
        REQUIRE(i->get_value() == 8);
        ++n_eights;
    }
    REQUIRE(n_eights == 2);

    unsigned n_nines = 0;
    for(auto i : c3.get_iter([] (int_t i) { return i==9;} )) {
        REQUIRE(i.get_ptr());
        REQUIRE(i->get_value() == 9);
        ++n_nines;
    }
    REQUIRE(n_nines == 2);

    unsigned n_nothing = 0;
    for(auto i : c3.get_iter([] (int_t) { return false;} )) {
        REQUIRE(i.get_ptr());
        ++n_nothing;
    }
    REQUIRE(n_nothing == 0);

    auto all_nines = c3.get_ptr_list([] (int_t i) { return i==9;});
    REQUIRE(all_nines.size() == 2);

    auto all_zeros = c3.get_ptr_list([] (int_t i) { return i;});
    REQUIRE(all_zeros.size() == 0);

    // number of emplace_back calls!
    REQUIRE(int_t::n_constructed == 3);

}

void TestRMSIQR() {

    // test exceptional cases
    {
        std_ext::IQR iqr;
        REQUIRE_THROWS_AS(iqr.GetMedian(), std::out_of_range);
        REQUIRE_THROWS_AS(iqr.GetIQRStdDev(), std::out_of_range);
    }

    // test simple, consecutive medians
    {
        std_ext::IQR iqr;
        iqr.Add(1);
        CHECK(iqr.GetMedian() == Approx(1));
        iqr.Add(2);
        iqr.Add(3);
        CHECK(iqr.GetMedian() == Approx(2));
        iqr.Add(4);
        CHECK(iqr.GetMedian() == Approx(2.5));
    }

    // compare to RMS of normal distribution
    {
        std::default_random_engine rnd;
        std::normal_distribution<double> gaussian(2,5);
        std_ext::RMS rms;
        std_ext::IQR iqr;
        for(auto i=0u;i<1e6;i++) {
            const auto v = gaussian(rnd);
            rms.Add(v);
            iqr.Add(v);
        }
        CHECK(rms.GetMean() == Approx(2).epsilon(0.01));
        CHECK(rms.GetRMS() == Approx(5).epsilon(0.01));
        CHECK(iqr.GetMedian() == Approx(rms.GetMean()).epsilon(0.01));
        CHECK(iqr.GetIQRStdDev() == Approx(rms.GetRMS()).epsilon(0.01));
        // check robustness by adding an outlier
        // RMS.GetMean changes alot, but IQR stays at correct value
        rms.Add(1e6);
        iqr.Add(1e6);
        CHECK(rms.GetMean()!=Approx(2).epsilon(0.01));
        CHECK(iqr.GetMedian()==Approx(2).epsilon(0.01));
    }
}


void TestDereference() {
    struct A {
        bool check() { return true; }
    };

    A a;
    auto a_shared = make_shared<A>(a);
    auto a_unique = std_ext::make_unique<A>(a);

    REQUIRE(std_ext::dereference(a).check());
    REQUIRE(std_ext::dereference(addressof(a)).check());
    REQUIRE(std_ext::dereference(a_shared).check());
    REQUIRE(std_ext::dereference(a_unique).check());
}
