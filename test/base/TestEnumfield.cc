#include "catch.hpp"

#include "base/enumfield.h"

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Enumfield", "[base]") {
    dotest();
}

enum class MyEnum {
    F1,F2
};

void dotest() {
    enumfield<MyEnum> field;

    REQUIRE_FALSE(field.isChecked(MyEnum::F1));
    REQUIRE_FALSE(field.isChecked(MyEnum::F2));

    REQUIRE_NOTHROW(field.Set(MyEnum::F2));
    REQUIRE(field.isChecked(MyEnum::F2));
    REQUIRE_FALSE(field.isChecked(MyEnum::F1));

    REQUIRE_NOTHROW(field.Reset(MyEnum::F2));
    REQUIRE_FALSE(field.isChecked(MyEnum::F2));
    REQUIRE_FALSE(field.isChecked(MyEnum::F1));

    REQUIRE_NOTHROW(field = MyEnum::F1);
    REQUIRE(field.isChecked(MyEnum::F1));
    REQUIRE(field == MyEnum::F1);
    const auto field2 = field;

    REQUIRE(field2==field);

    constexpr enumfield<MyEnum> cf1(MyEnum::F1);
    constexpr enumfield<MyEnum> cf2(MyEnum::F2);

    REQUIRE(cf1.isChecked(MyEnum::F1));

    constexpr auto cf3 = cf1 | cf2;
    REQUIRE(cf3.isChecked(MyEnum::F1));
    REQUIRE(cf3.isChecked(MyEnum::F2));

    static_assert(cf3.isChecked(MyEnum::F1),"constexpr isChecked");

    static_assert(std::is_default_constructible<enumfield<MyEnum>>::value,"");
    static_assert(std::is_move_assignable<enumfield<MyEnum>>::value,"");
    static_assert(std::is_move_constructible<enumfield<MyEnum>>::value,"");



}
