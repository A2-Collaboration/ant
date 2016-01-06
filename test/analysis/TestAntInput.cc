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

    TEvent* event = new TEvent();

    //======= Candidate 1 =============

    event->Candidates.push_back(
          TCandidate(Detector_t::Type_t::CB | Detector_t::Type_t::PID, 200,2,2,1,2,0,0,{})
          );

    event->Candidates.back().Clusters.push_back(
          TCluster(TVector3(25,0,0), 270, -260, Detector_t::Type_t::CB,0)
          );


    event->Candidates.back().Clusters.back().Hits.push_back(
          TClusterHit(110, {
                             {Channel_t::Type_t::Integral, 150}, // MeV
                             {Channel_t::Type_t::Timing, -290}   // ns
                           })
          );

    event->Candidates.back().Clusters.back().Hits.push_back(
          TClusterHit(220, {
                             {Channel_t::Type_t::Integral, 120}, // MeV
                             {Channel_t::Type_t::Timing, -280}   // ns
                           })
          );

    event->Candidates.back().Clusters.push_back(
          TCluster(TVector3(10,0,0), 5, -260, Detector_t::Type_t::PID,0)
          );

    event->Candidates.back().Clusters.back().Hits.push_back(
          TClusterHit(20, {{Channel_t::Type_t::Integral, 4}})
          );

    //======= Candidate 2 =============
    event->Candidates.push_back(
                TCandidate(Detector_t::Type_t::TAPS | Detector_t::Type_t::TAPSVeto, 300,2,3,0.5,0,0,0,{})
          );
    event->Candidates.back().Clusters.push_back(
          TCluster(TVector3(25,10,0),300, -160, Detector_t::Type_t::TAPS,0)
          );
    event->Candidates.back().Clusters.push_back(
          TCluster(TVector3(10,1,0),5,-150, Detector_t::Type_t::TAPSVeto,0)
          );

    cout << *event << endl;

    auto antevent = analysis::input::Converter::Convert(*event);

    cout << antevent << endl;

    const analysis::data::CandidateList& cand = antevent.Reconstructed.Candidates;

    REQUIRE(cand.size() == 2);
    /// \todo Write more REQUIRE stuff here

    REQUIRE(cand.at(0)->ClusterSize==2);
    REQUIRE(cand.at(1)->ClusterSize==0);

    REQUIRE(cand.at(0)->GetDetector()==(Detector_t::Type_t::CB | Detector_t::Type_t::PID));
    REQUIRE(cand.at(1)->GetDetector()==(Detector_t::Type_t::TAPS | Detector_t::Type_t::TAPSVeto));

}


