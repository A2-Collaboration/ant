#include "catch.hpp"

#include "base/Detector_t.h"

#include <sstream>

using namespace ant;
using namespace std;

const Detector_t::Any_t CBTAPS = Detector_t::Type_t::CB | Detector_t::Type_t::TAPS;

Detector_t::Any_t cb_taps_test(const Detector_t::Any_t& a, const Detector_t::Any_t& b) {
    return (a & CBTAPS) ^ (b & CBTAPS);
}

TEST_CASE("Detector_t: Any_t", "[base]") {

    // test printable
    stringstream ss;
    ss << Detector_t::Any_t::CB_Apparatus;
    REQUIRE(ss.str() == "CB|PID|MWPC0|MWPC1");

    // test easy usage in if-statements
    if(Detector_t::Any_t::CB_Apparatus & Detector_t::Type_t::CB) {
        REQUIRE(true);
    }
    else {
        REQUIRE(false);
    }

    REQUIRE((Detector_t::Any_t::CB_Apparatus & Detector_t::Type_t::CB));
    REQUIRE((Detector_t::Type_t::MWPC0 & Detector_t::Any_t::CB_Apparatus));
    REQUIRE((Detector_t::Type_t::MWPC1 & Detector_t::Any_t::CB_Apparatus));
    REQUIRE((Detector_t::Type_t::PID & Detector_t::Any_t::CB_Apparatus));
    REQUIRE_FALSE((Detector_t::Type_t::TAPS & Detector_t::Any_t::CB_Apparatus));
    REQUIRE_FALSE((Detector_t::Type_t::TAPSVeto & Detector_t::Any_t::CB_Apparatus));

    REQUIRE((Detector_t::Type_t::TAPS & Detector_t::Any_t::TAPS_Apparatus));
    REQUIRE((Detector_t::Type_t::TAPSVeto & Detector_t::Any_t::TAPS_Apparatus));

    REQUIRE((Detector_t::Type_t::PID & Detector_t::Any_t::Veto));
    REQUIRE((Detector_t::Type_t::TAPSVeto & Detector_t::Any_t::Veto));

    REQUIRE(Detector_t::Any_t::CB_Apparatus == Detector_t::Any_t::CB_Apparatus);
    REQUIRE(Detector_t::Any_t::CB_Apparatus != Detector_t::Any_t::TAPS_Apparatus);

    Detector_t::Any_t detector1 = Detector_t::Any_t::None;
    detector1 |= Detector_t::Type_t::CB;

    REQUIRE((detector1 & Detector_t::Any_t::CB_Apparatus));
    REQUIRE(detector1 == Detector_t::Type_t::CB);
    REQUIRE(detector1 != Detector_t::Type_t::TAPS);
    REQUIRE(detector1 != Detector_t::Type_t::MWPC0);

    REQUIRE((cb_taps_test(Detector_t::Type_t::CB, Detector_t::Type_t::TAPS) & CBTAPS));
    REQUIRE((cb_taps_test(Detector_t::Type_t::TAPS, Detector_t::Type_t::CB) & CBTAPS));
    REQUIRE_FALSE((cb_taps_test(Detector_t::Type_t::CB, Detector_t::Type_t::CB) & CBTAPS));
    REQUIRE_FALSE((cb_taps_test(Detector_t::Type_t::TAPS, Detector_t::Type_t::TAPS) & CBTAPS));

    REQUIRE((cb_taps_test(Detector_t::Type_t::CB | Detector_t::Type_t::PID, Detector_t::Type_t::TAPS) & CBTAPS));

    REQUIRE(Detector_t::Any_t::CB_Apparatus.test(Detector_t::Type_t::CB));
}

TEST_CASE("Detector_t: ElementFlags_t", "[base]") {

    Detector_t::ElementFlags_t flags;
    flags.set(Detector_t::ElementFlag_t::BadTDC);
    flags.set(Detector_t::ElementFlag_t::Missing);
    REQUIRE(flags.test(Detector_t::ElementFlag_t::Missing));
    REQUIRE(flags.test(Detector_t::ElementFlag_t::BadTDC));
    flags.unset(Detector_t::ElementFlag_t::Missing);
    REQUIRE_FALSE(flags.test(Detector_t::ElementFlag_t::Missing));
    REQUIRE(flags.test(Detector_t::ElementFlag_t::BadTDC));

    REQUIRE(flags.any());
    REQUIRE_FALSE(flags.none());
    REQUIRE_FALSE(flags.all());
    REQUIRE(static_cast<bool>(flags));
}

struct TestDetector_t : Detector_t {
    std::vector<Detector_t::Element_t> Elements;

    TestDetector_t() : Detector_t(Detector_t::Type_t::CB) {
        Elements.resize(10, Element_t(0, {0,0,0}));
    }

    virtual unsigned GetNChannels() const override
    {
        return Elements.size();
    }
    virtual vec3 GetPosition(unsigned channel) const override
    {
        return Elements[channel].Position;
    }
    virtual void SetElementFlags(unsigned channel, const ElementFlags_t& flags) override
    {
        Elements[channel].Flags |= flags;
    }
    virtual const ElementFlags_t& GetElementFlags(unsigned channel) const override
    {
        return Elements[channel].Flags;
    }
};

TEST_CASE("Detector_t: Detector interface", "[base]") {
    TestDetector_t det;
    CHECK(det.Type == Detector_t::Type_t::CB);
    CHECK(static_cast<bool>(det.Type & Detector_t::Any_t::Calo));
    CHECK(static_cast<bool>(det.Type & Detector_t::Any_t::CB_Apparatus));

    CHECK(!det.IsIgnored(0));
    det.SetElementFlags(0, Detector_t::ElementFlag_t::Broken);
    CHECK(det.IsIgnored(0));

    CHECK(!det.IsIgnored(1));
    det.SetElementFlags(1, Detector_t::ElementFlag_t::Missing);
    CHECK(det.IsIgnored(1));

    CHECK(!det.IsIgnored(2));
    det.SetElementFlags(2, Detector_t::ElementFlag_t::Missing);
    det.SetElementFlags(2, Detector_t::ElementFlag_t::Broken);
    CHECK(det.IsIgnored(2));
    CHECK(!det.HasElementFlags(2, Detector_t::ElementFlag_t::BadTDC));

    det.SetElementFlag(Detector_t::ElementFlag_t::Broken, {4, 6, 8});
    CHECK(det.IsIgnored(4));
    CHECK(det.IsIgnored(6));
    CHECK(det.IsIgnored(8));
}

