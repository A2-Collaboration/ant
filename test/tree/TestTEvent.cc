#include "catch.hpp"
#include "tree/TEvent.h"
#include "tree/TTrack.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "TRandom2.h"

using namespace std;
using namespace ant;

TRandom2 rng;


TEST_CASE("TEvent: Write to TTree", "[tree]") {

    TFile f("TEvent_test.root","RECREATE");

    TTree* tree = new TTree("teventtest","TEvent Test Tree");
    TEvent* buffer = new TEvent();

    tree->Branch("event", buffer);

    //======= Track 1 =============

    LogicalChannel_t lc;


    buffer->Tracks.push_back(
                TTrack(200,1,2,2,0,0)
            );

    buffer->Tracks.back().Clusters.push_back(
                TCluster(TVector3(25,0,0),200,Detector_t::Type_t::CB)
                );

    lc.Channel = 1;
    lc.ChannelType = Channel_t::Type_t::Integral;
    lc.DetectorType = Detector_t::Type_t::CB;

    buffer->Tracks.back().Clusters.back().Hits.push_back(
                TDetectorReadHit(lc,{2})
                );

    lc.Channel = 2;
    lc.ChannelType = Channel_t::Type_t::Integral;
    lc.DetectorType = Detector_t::Type_t::CB;
    buffer->Tracks.back().Clusters.back().Hits.push_back(
                TDetectorReadHit(lc,{3})
                );

    buffer->Tracks.back().Clusters.push_back(
                TCluster(TVector3(10,0,0),5,Detector_t::Type_t::PID)
                );

    lc.Channel = 20;
    lc.ChannelType = Channel_t::Type_t::Integral;
    lc.DetectorType = Detector_t::Type_t::PID;
    buffer->Tracks.back().Clusters.back().Hits.push_back(
                TDetectorReadHit(lc,{6})
                );
    lc.Channel = 20;
    lc.ChannelType = Channel_t::Type_t::Timing;
    lc.DetectorType = Detector_t::Type_t::PID;
    buffer->Tracks.back().Clusters.back().Hits.push_back(
                TDetectorReadHit(lc,{1})
                );

    //======= Track 1 =============
    buffer->Tracks.push_back(
                TTrack(300,0.5,2,3,0,0)
            );
    buffer->Tracks.back().Clusters.push_back(
                TCluster(TVector3(25,10,0),300,Detector_t::Type_t::TAPS)
                );
    buffer->Tracks.back().Clusters.push_back(
                TCluster(TVector3(10,1,0),5,Detector_t::Type_t::TAPSVeto)
                );

    tree->Fill();

    cout << buffer << endl;

    f.Write();
    f.Close();

    tree = nullptr;

    TFile f2("TEvent_test.root","READ");
    REQUIRE(f2.IsOpen());

    f2.GetObject("teventtest", tree);

    REQUIRE(tree!=nullptr);

    TEvent* readback = nullptr;

    tree->SetBranchAddress("event",&readback);

    REQUIRE(tree->GetEntries()==1);

    tree->GetEntry(0);

    cout << *readback << endl;

    //TODO: Compate the events. implement == operator for events...


}
