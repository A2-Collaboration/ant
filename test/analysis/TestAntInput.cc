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

    //======= Candidate 1 =============

    event->Candidates.push_back(
          ant::TCandidate(200,1,2,2,{})
          );

    event->Candidates.back().Clusters.push_back(
          ant::TCluster(TVector3(25,0,0), 270, -260, 7, ant::Detector_t::Type_t::CB,0)
          );


    event->Candidates.back().Clusters.back().Hits.push_back(
          ant::TClusterHit(110, {
                             {ant::Channel_t::Type_t::Integral, 150}, // MeV
                             {ant::Channel_t::Type_t::Timing, -290}   // ns
                           })
          );

    event->Candidates.back().Clusters.back().Hits.push_back(
          ant::TClusterHit(220, {
                             {ant::Channel_t::Type_t::Integral, 120}, // MeV
                             {ant::Channel_t::Type_t::Timing, -280}   // ns
                           })
          );

    event->Candidates.back().Clusters.push_back(
          ant::TCluster(TVector3(10,0,0), 5, -260, 7, ant::Detector_t::Type_t::PID,0)
          );

    event->Candidates.back().Clusters.back().Hits.push_back(
          ant::TClusterHit(20, {{ant::Channel_t::Type_t::Integral, 4}})
          );

    //======= Candidate 2 =============
    event->Candidates.push_back(
          ant::TCandidate(300,0.5,2,3)
          );
    event->Candidates.back().Clusters.push_back(
          ant::TCluster(TVector3(25,10,0),300, -160,7, ant::Detector_t::Type_t::TAPS,0)
          );
    event->Candidates.back().Clusters.push_back(
          ant::TCluster(TVector3(10,1,0),5,-150,9, ant::Detector_t::Type_t::TAPSVeto,0)
          );

    cout << *event << endl;

    auto antevent = ant::input::Converter::Convert(*event);

    cout << antevent << endl;

    const CandidateList& cand = antevent.Reconstructed().Candidates();

    REQUIRE(cand.size() == 2);
    /// \todo Write more REQUIRE stuff here

    REQUIRE(cand.at(0)->ClusterSize()==2);
    REQUIRE(cand.at(1)->ClusterSize()==0);

    REQUIRE(cand.at(0)->Detector()==(detector_t::CB | detector_t::PID));
    REQUIRE(cand.at(1)->Detector()==(detector_t::TAPS | detector_t::TAPSVeto));

}


