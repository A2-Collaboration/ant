#include "catch.hpp"

#include "tree/TEvent.h"
#include "tree/TCandidate.h"

#include "base/tmpfile_t.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

void dotest();

TEST_CASE("TEvent: Write to TTree", "[tree]") {
  dotest();
}

void dotest() {
  ant::tmpfile_t tmpfile;


  TFile f(tmpfile.filename.c_str(),"RECREATE");

  TTree* tree = new TTree("teventtest","TEvent Test Tree");
  ant::TEvent* buffer = new ant::TEvent();

  tree->Branch("event", buffer);

  //======= Track 1 =============



  buffer->Candidates.push_back(
        ant::TCandidate(200,1,2,2,{})
        );

  buffer->Candidates.back().Clusters.push_back(
        ant::TCluster(TVector3(25,0,0), 270, -350, ant::Detector_t::Type_t::CB)
        );


  buffer->Candidates.back().Clusters.back().Hits.push_back(
        ant::TClusterHit(110, {
                           {ant::Channel_t::Type_t::Integral, 150}, // MeV
                           {ant::Channel_t::Type_t::Timing, -290}   // ns
                         })
        );

  buffer->Candidates.back().Clusters.back().Hits.push_back(
        ant::TClusterHit(220, {
                           {ant::Channel_t::Type_t::Integral, 120}, // MeV
                           {ant::Channel_t::Type_t::Timing, -280}   // ns
                         })
        );

  buffer->Candidates.back().Clusters.push_back(
        ant::TCluster(TVector3(10,0,0), 5, -270, ant::Detector_t::Type_t::PID)
        );

  buffer->Candidates.back().Clusters.back().Hits.push_back(
        ant::TClusterHit(20, {{ant::Channel_t::Type_t::Integral, 4}})
        );

  //======= Track 2 =============
  buffer->Candidates.push_back(
        ant::TCandidate(300,0.5,2,3)
        );
  buffer->Candidates.back().Clusters.push_back(
        ant::TCluster(TVector3(25,10,0),300,-200,ant::Detector_t::Type_t::TAPS)
        );
  buffer->Candidates.back().Clusters.push_back(
        ant::TCluster(TVector3(10,1,0),5,-250, ant::Detector_t::Type_t::TAPSVeto)
        );

  tree->Fill();

  cout << buffer << endl;

  f.Write();
  f.Close();

  tree = nullptr;

  TFile f2(tmpfile.filename.c_str(),"READ");
  REQUIRE(f2.IsOpen());

  f2.GetObject("teventtest", tree);

  REQUIRE(tree!=nullptr);

  ant::TEvent* readback = nullptr;

  tree->SetBranchAddress("event",&readback);

  REQUIRE(tree->GetEntries()==1);

  tree->GetEntry(0);

  cout << *readback << endl;

  /// \todo Compare the events. implement == operator for events...


}
