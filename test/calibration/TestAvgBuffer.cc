#include "catch.hpp"

#include "calibration/gui/AvgBuffer.h"

#include "tree/TID.h"

#include "base/std_ext/string.h"
#include "base/interval.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::calibration::gui;

template<typename Hist = TH1D>
shared_ptr<Hist> makeHist(double value) {
    static unsigned num = 0;
    string title = std_ext::formatter() << "hist_" << num++;
    auto hist = std::make_shared<Hist>(title.c_str(), title.c_str(), 1, 0, 1);
    hist->SetDirectory(0);
    hist->Fill(0.0, value);
    return hist;
}

TEST_CASE("TestAvgBuffer: AvgBuffer_Sum","[calibration]")
{
    AvgBuffer_Sum buf;
    buf.Push(makeHist(1), {1,1});
    buf.Push(makeHist(2), {2,2});
    buf.Push(makeHist(3), {3,3});
    buf.Push(makeHist(4), {4,4});
    // stays empty no matter how many histograms you push
    // until calling Flush()
    REQUIRE(buf.Empty());
    buf.Flush();
    REQUIRE_FALSE(buf.Empty());
    REQUIRE(buf.CurrentHist().GetBinContent(1) == 10);
    // becomes empty after first Next()
    buf.Next();
    REQUIRE(buf.Empty());
}

void dotest_savitzkygolay_simple();
void dotest_savitzkygolay_avg();

TEST_CASE("TestAvgBuffer: AvgBuffer_SavitzkyGolay simple","[calibration]") {
    dotest_savitzkygolay_simple();
}

TEST_CASE("TestAvgBuffer: AvgBuffer_SavitzkyGolay average","[calibration]") {
    dotest_savitzkygolay_avg();
}

void dotest_savitzkygolay_simple()
{
    AvgBuffer_SavitzkyGolay buf(5,4);
    for(int i=0;i<20;i++)
        buf.Push(makeHist(1), {uint32_t(i),uint32_t(i)});
    buf.Flush();
    unsigned nNext = 0;
    while(!buf.Empty()) {
        INFO(nNext++);
        REQUIRE(buf.CurrentHist().GetBinContent(1)==Approx(1));
        buf.Next();
    }
    REQUIRE(nNext==20);
}

vector<double> calc_moving_avg(const vector<double>& data, unsigned avgLength)
{
    if(avgLength<1)
        avgLength = 1;

    const unsigned n = data.size();
    if(avgLength>n)
        avgLength=n;

    auto get_data = [&data] (const int i) {
        if(i<0)
            return data.at(-i);
        if(unsigned(i)>=data.size()) {
            const auto i_ = data.size()-(i-data.size())-2;
            return data.at(i_);
        }
        return data[i];
    };

    const int n_l = (int(avgLength)-1)/2 + (avgLength % 2 == 0);
    const int n_r = (int(avgLength)-1)/2;

    vector<double> moving_avg(n, 0); // init with zeros
    for(unsigned i=0;i<n;i++) {

        for(int j=-n_l;j<=n_r;j++) {
            moving_avg[i] += get_data(i+j);
        }
        moving_avg[i] /= double(avgLength);
    }
    return moving_avg;
}

void dotest_savitzkygolay_avg()
{
    const vector<double> data = {1, 32, 12, 46, 15, 61, 3, 4, 10, 13,
                                 4, 5, 6, 9, 1, 10, 11, 18, 39, 10};

    const string data_str = std_ext::formatter() << data;
    INFO("Data=" << data_str);

    for(unsigned avgLength=1;avgLength<=data.size();avgLength++) {
        INFO("avgLength=" << avgLength);

        AvgBuffer_SavitzkyGolay buf(avgLength, 0); // order 0 should be moving average
        unsigned nNextID = 0;
        unsigned nPushed = 0;
        const vector<double> expected = calc_moving_avg(data, avgLength);
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
                auto value = buf.CurrentHist().GetBinContent(1);
                REQUIRE(value == Approx(expected[nNextID]));
                REQUIRE(buf.CurrentRange() == interval<TID>(nNextID, nNextID));
                nNextID++;
                buf.Next();
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
