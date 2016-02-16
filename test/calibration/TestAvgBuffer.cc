#include "catch.hpp"

#include "calibration/gui/AvgBuffer.h"

#include "tree/TID.h"

#include "base/std_ext/string.h"
#include "base/interval.h"

#include "TH1D.h"

#include <iostream>

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
    buf.Push(makeHist(2), {2,2});
    buf.Push(makeHist(3), {3,3});
    buf.Push(makeHist(4), {4,4});
    // avgLength==0 is always empty until flush
    REQUIRE(buf.Empty());
    buf.Flush();
    REQUIRE_FALSE(buf.Empty());
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 10);
}


TEST_CASE("TestAvgBuffer: AvgLength=1","[calibration]")
{
    AvgBuffer<TH1D, interval<int> > buf(1);

    buf.Push(makeHist(1), {1,1});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 1);
    buf.GotoNextID();

    buf.Push(makeHist(2), {2,2});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 2);
    buf.GotoNextID();

    buf.Push(makeHist(3), {3,3});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 3);
    buf.GotoNextID();

    buf.Push(makeHist(4), {4,4});
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 4);
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
    buf.GotoNextID();

    buf.Push(makeHist(3), {3,3});
    buf.Push(makeHist(4), {4,4});
    buf.GotoNextID();
    buf.GotoNextID();
    REQUIRE_FALSE(buf.Empty());
    REQUIRE(buf.CurrentSum()->GetBinContent(1) == 7);
    buf.GotoNextID();
    REQUIRE(buf.Empty());
}

void dotest();

TEST_CASE("TestAvgBuffer: Complete run","[calibration]")
{
    dotest();
}

vector<double> calc_moving_sum(const vector<double>& data, unsigned avgLength)
{
    if(avgLength<1)
        avgLength = 1;

    const unsigned n = data.size();
    if(avgLength>n)
        avgLength=n;

    vector<double> moving_sum(n, 0); // init with zeros
    for(unsigned i=0;i<n;i++) {
        interval<unsigned> window{i-avgLength/2, i + avgLength/2 + avgLength%2};

        if(i<=avgLength/2)
            window = {0, avgLength};
        if(window.Stop() > n)
            window = {n-avgLength, n};

        for(unsigned j=window.Start();j<window.Stop();j++)
            moving_sum[i] += data[j];
    }
    return moving_sum;
}

void dotest()
{
    const vector<double> data = {1, 32, 12, 46, 15, 61, 3, 4, 10, 13,
                                 4, 5, 6, 9, 1, 10, 11, 18, 39, 10};

    const string data_str = std_ext::formatter() << data;
    INFO("Data=" << data_str);

    for(unsigned avgLength=1;avgLength<=data.size();avgLength++) {
        INFO("avgLength=" << avgLength);

        AvgBuffer<TH1D, interval<unsigned> > buf(avgLength);
        unsigned nNextID = 0;
        unsigned nPushed = 0;
        const vector<double> expected = calc_moving_sum(data, avgLength);
        const string expected_str = std_ext::formatter() << expected;
        INFO("Expected=" << expected_str);

        list<pair<unsigned, double>> buf_items;
        for(unsigned i=0;i<data.size();i++)
            buf_items.emplace_back(i, data[i]);

        do {

            while(buf.Empty() && !buf_items.empty()) {
                const unsigned id = buf_items.front().first;
                const double value = buf_items.front().second;
                buf.Push(makeHist(value),{id,id});
                buf_items.pop_front();
                nPushed++;
            }

            INFO("nPushed=" << nPushed);

            if(!buf.Empty()) {
                INFO("nNextID=" << nNextID);
                auto value = buf.CurrentSum()->GetBinContent(1);
                REQUIRE(value == expected[nNextID]);
                REQUIRE(buf.CurrentID() == interval<unsigned>(nNextID, nNextID));
                nNextID++;
                buf.GotoNextID();
            }
            if(buf_items.empty()) {
                buf.Flush();
            }
        }
        while(!buf_items.empty() || !buf.Empty());

        REQUIRE(nNextID == data.size());
        REQUIRE(nPushed == data.size());

    }
}

