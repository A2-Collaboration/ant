#include "catch.hpp"

#include "base/OptionsList.h"

#include <cassert>
#include <iostream>


using namespace std;
using namespace ant;

void dotest1();
void dotest2();
void dotest3();
void dotest4();
void dotest5();

TEST_CASE("OptionsList: Basic", "[base]") {
    dotest1();
}

TEST_CASE("OptionsList: Chained", "[base]") {
    dotest2();
}

TEST_CASE("OptionsList: Flags", "[base]") {
    dotest3();
}

TEST_CASE("OptionsList: Unused/Notfound", "[base]") {
    dotest4();
}

TEST_CASE("OptionsList: Has (unused) strings", "[base]") {
    dotest5();
}


void dotest1() {
    std::shared_ptr<OptionsList> opts;

    REQUIRE_NOTHROW(opts = std::make_shared<OptionsList>());
    REQUIRE_NOTHROW(opts->SetOption("key=val"));

    REQUIRE(opts->Get<string>("key")  == "val");
    REQUIRE(opts->Get<string>("key2") == "");
}

void dotest2() {
    auto opts = std::make_shared<OptionsList>();
    opts->SetOption("key=val");

    std::shared_ptr<OptionsList> opts2;
    REQUIRE_NOTHROW(opts2 = std::make_shared<OptionsList>(opts));

    REQUIRE_NOTHROW(opts2->SetOption("key4=val4"));    // own option
    REQUIRE(opts2->Get<string>("key4") == "val4");

    REQUIRE(opts2->Get<string>("key")  == "val"); // from parent
    REQUIRE(opts2->Get<string>("key2") == "");    // from parent

    REQUIRE_NOTHROW(opts2->SetOption("key2=val2"));     // overwrite parent
    REQUIRE(opts2->Get<string>("key2") == "val2");
}

void dotest3() {
    auto opts = std::make_shared<OptionsList>();
    opts->SetOption("flag1=1");
    opts->SetOption("flag2=on");
    opts->SetOption("flag3=oN");
    opts->SetOption("flag4=truE");
    opts->SetOption("flag5=no");
    opts->SetOption("flag6=YES");
    REQUIRE(opts->Get<bool>("flag1"));
    REQUIRE(opts->Get<bool>("flag2"));
    REQUIRE(opts->Get<bool>("flag3"));
    REQUIRE(opts->Get<bool>("flag4"));
    REQUIRE_FALSE(opts->Get<bool>("flag5"));
    REQUIRE(opts->Get<bool>("flag6"));
}

void dotest4() {
    auto opts = std::make_shared<OptionsList>();
    opts->SetOption("flag1=1");
    opts->SetOption("flag2=on");
    opts->SetOption("flag3=oN");
    opts->SetOption("flag4=truE");
    opts->SetOption("flag5=no");
    opts->SetOption("flag6=YES");

    REQUIRE(opts->Get<bool>("flag1"));
    REQUIRE(opts->Get<bool>("flag3"));
    REQUIRE(opts->Get<bool>("flag7", true));

    auto notfound = opts->GetNotFound();
    REQUIRE(notfound.size() == 1);

    auto unused = opts->GetUnused();
    REQUIRE(unused.size()==4);
}

void dotest5() {
    auto opts = std::make_shared<OptionsList>();
    opts->SetOption("someOptionString=value");
    opts->SetOption("some_other_string=3");
    opts->SetOption("moar_strings=");

    std::shared_ptr<OptionsList> child_opts;
    REQUIRE_NOTHROW(child_opts = std::make_shared<OptionsList>(opts));
    child_opts->SetOption("key=val");

    REQUIRE(opts->HasOptionStartsWith("someOption"));
    REQUIRE_FALSE(opts->HasOptionStartsWith("string"));
    string key;
    REQUIRE_NOTHROW(key = child_opts->OptionStartsWith("someOption"));
    REQUIRE(child_opts->Get<string>(key) == "value");
    REQUIRE(opts->Get<double>(opts->OptionStartsWith("some_other")) == 3);
    REQUIRE_FALSE(opts->HasUnusedOptionStartsWith("some"));
    REQUIRE(opts->HasUnusedOptionStartsWith("moar"));
}
