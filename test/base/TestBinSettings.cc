#include "catch.hpp"
#include "base/BinSettings.h"
#include "base/std_ext/math.h"
#include "base/std_ext/memory.h"
#include <iostream>
#include <sstream>

#include "tclap/CmdLine.h"


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
    ss_b << "(200,[0.95:1.2])";
    REQUIRE(ss_b >> b);
    REQUIRE(b.Bins()==200);
    REQUIRE(b.Start()==0.95);
    REQUIRE(b.Stop()==1.2);
}

TEST_CASE("BinSettings: Parse from istringstream (TCLAP test code)", "[base]")
{
    BinSettings destVal(0);

    const string strVal("(200,[0.95:1.2])");
    std::istringstream is(strVal);

    int valuesRead = 0;
    while ( is.good() ) {
        if ( is.peek() != EOF )
            is >> destVal;
        else
            break;
        valuesRead++;
    }
    REQUIRE_FALSE(is.fail());
}

struct TCLAPBinSettings : BinSettings {
    using BinSettings::BinSettings;
    using ValueCategory = TCLAP::ValueLike;
};


TEST_CASE("BinSettings: Use as cmdline parameter", "[base]") {
    TCLAP::CmdLine cmd("TestCmdLine");
    cmd.setExceptionHandling(false);
    auto cmd_test = cmd.add<TCLAP::ValueArg<TCLAPBinSettings>>("" ,"test", "", false, TCLAPBinSettings(0), "");

    vector<string> args{"progname", "--test", "(200,[0.95:1.2])"};
    REQUIRE_NOTHROW(cmd.parse(args));
}
