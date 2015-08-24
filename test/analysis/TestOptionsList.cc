#include "catch.hpp"

#include "analysis/physics/Physics.h"

#include <cassert>
#include <iostream>


using namespace std;
using namespace ant;
using namespace ant::analysis;

void dotest1();
void dotest2();

TEST_CASE("optionsList: Basic", "[analysis]") {
    dotest1();
}

TEST_CASE("optionsList: Chained", "[analysis]") {
    dotest2();
}
std::shared_ptr<OptionsList> opts;

void dotest1() {

    REQUIRE_NOTHROW(opts = std::make_shared<OptionsList>());
    REQUIRE_NOTHROW(opts->SetOption("key=val"));

    REQUIRE(opts->GetOption("key")  == string("val"));
    REQUIRE(opts->GetOption("key2") == string(""));
}

void dotest2() {
    std::shared_ptr<OptionsList> opts2;
    REQUIRE_NOTHROW(opts2 = std::make_shared<OptionsList>(opts));

    REQUIRE_NOTHROW(opts2->SetOption("key4=val4"));    // own option
    REQUIRE(opts2->GetOption("key4") == string("val4"));

    REQUIRE(opts2->GetOption("key")  == string("val")); // from parent
    REQUIRE(opts2->GetOption("key2") == string(""));    // from parent

    REQUIRE_NOTHROW(opts2->SetOption("key2=val2"));     // overwrite parent
    REQUIRE(opts2->GetOption("key2") == string("val2"));

}
