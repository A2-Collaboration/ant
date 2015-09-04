#include "catch.hpp"

#include "TCluster.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("TCluster","[tree]")
{
    dotest();
}

void dotest()
{
    TCluster cl;
    REQUIRE_FALSE(cl.HasFlag(TCluster::Flags_t::Split));
    REQUIRE_FALSE(cl.HasFlag(TCluster::Flags_t::TouchesHole));

    cl.SetFlag(TCluster::Flags_t::Split);
    REQUIRE(cl.HasFlag(TCluster::Flags_t::Split));
    REQUIRE_FALSE(cl.HasFlag(TCluster::Flags_t::TouchesHole));

    cl.SetFlag(TCluster::Flags_t::TouchesHole);
    REQUIRE(cl.HasFlag(TCluster::Flags_t::Split));
    REQUIRE(cl.HasFlag(TCluster::Flags_t::TouchesHole));

    cl.SetFlag(TCluster::Flags_t::Split, false);
    REQUIRE_FALSE(cl.HasFlag(TCluster::Flags_t::Split));
    REQUIRE(cl.HasFlag(TCluster::Flags_t::TouchesHole));
}
