#include "catch.hpp"

#include "calibration/gui/AvgBuffer.h"

#include "tree/TID.h"

#include "base/std_ext/string.h"
#include "base/interval.h"

#include "TH1D.h"

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;

shared_ptr<TH1D> makeHist(double value) {
    static unsigned num = 0;
    string title = std_ext::formatter() << "hist_" << num++;
    auto hist = std::make_shared<TH1D>(title.c_str(), title.c_str(), 1, 0, 1);
    hist->SetDirectory(0);
    hist->Fill(0.0, value);
    return hist;
}

TEST_CASE("TestAvgBuffer: AvgLength=0","[calibration]")
{
    AvgBuffer<TH1D, interval<int> > buf(0);
    buf.Push(makeHist(1), {1,1});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 1);
    buf.Push(makeHist(2), {2,2});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 3);
    buf.Push(makeHist(3), {3,3});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 6);
    buf.Push(makeHist(4), {4,4});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 10);
    // avgLength==0 is always empty
    REQUIRE(buf.Empty());
}


TEST_CASE("TestAvgBuffer: AvgLength=1","[calibration]")
{
    AvgBuffer<TH1D, interval<int> > buf(1);
    buf.Push(makeHist(1), {1,1});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 1);
    buf.Push(makeHist(2), {2,2});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 2);
    buf.Push(makeHist(3), {3,3});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 3);
    buf.Push(makeHist(4), {4,4});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 4);
    buf.GotoNextID();
    buf.GotoNextID();
    buf.GotoNextID();
    REQUIRE(buf.Empty());
}

TEST_CASE("TestAvgBuffer: AvgLength=2","[calibration]")
{
    AvgBuffer<TH1D, interval<int> > buf(2);
    buf.Push(makeHist(1), {1,1});
    REQUIRE(buf.Empty());
    buf.Push(makeHist(2), {2,2});
    REQUIRE_FALSE(buf.Empty());
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 3);

    buf.Push(makeHist(3), {3,3});
    buf.Push(makeHist(4), {4,4});
    REQUIRE_FALSE(buf.Empty());
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 7);
    buf.GotoNextID();
    buf.GotoNextID();
    buf.GotoNextID();
    REQUIRE(buf.Empty());
}




