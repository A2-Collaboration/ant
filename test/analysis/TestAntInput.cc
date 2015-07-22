#include "catch.hpp"

#include "analysis/input/ant/detail/Convert.h"
#include "analysis/data/Event.h"

#include "tree/TEvent.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("AntInput", "[analysis]") {
    dotest();
}

void dotest() {

    ant::TEvent* event = new ant::TEvent();

    //======= Track 1 =============

    event->Tracks.push_back(
          ant::TTrack(200,1,2,2,{})
          );

    event->Tracks.back().Clusters.push_back(
          ant::TCluster(TVector3(25,0,0), 270, ant::Detector_t::Type_t::CB)
          );


    event->Tracks.back().Clusters.back().Hits.push_back(
          ant::TClusterHit(110, {
                             {ant::Channel_t::Type_t::Integral, 150}, // MeV
                             {ant::Channel_t::Type_t::Timing, -290}   // ns
                           })
          );

    event->Tracks.back().Clusters.back().Hits.push_back(
          ant::TClusterHit(220, {
                             {ant::Channel_t::Type_t::Integral, 120}, // MeV
                             {ant::Channel_t::Type_t::Timing, -280}   // ns
                           })
          );

    event->Tracks.back().Clusters.push_back(
          ant::TCluster(TVector3(10,0,0), 5, ant::Detector_t::Type_t::PID)
          );

    event->Tracks.back().Clusters.back().Hits.push_back(
          ant::TClusterHit(20, {{ant::Channel_t::Type_t::Integral, 4}})
          );

    //======= Track 2 =============
    event->Tracks.push_back(
          ant::TTrack(300,0.5,2,3)
          );
    event->Tracks.back().Clusters.push_back(
          ant::TCluster(TVector3(25,10,0),300,ant::Detector_t::Type_t::TAPS)
          );
    event->Tracks.back().Clusters.push_back(
          ant::TCluster(TVector3(10,1,0),5, ant::Detector_t::Type_t::TAPSVeto)
          );

    cout << *event << endl;

    shared_ptr<ant::Event> antevent = ant::input::Convert(*event);

    cout << *antevent << endl;

    REQUIRE(antevent->Reconstructed().Tracks().size() == 2);
    /// \todo Write more REQUIRE stuff here

}


