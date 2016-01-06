#include "catch.hpp"

#include "tree/TEvent.h"
#include "tree/TCandidate.h"

#include "base/tmpfile_t.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("TEvent: Write to TTree", "[tree]") {
  dotest();
}

void dotest() {
  tmpfile_t tmpfile;


  TFile f(tmpfile.filename.c_str(),"RECREATE");

  TTree* tree = new TTree("teventtest","TEvent Test Tree");
  TEvent* buffer = new TEvent();

  tree->Branch("event", buffer);

  //======= Candidate 1 =============



  buffer->Candidates.push_back(
        TCandidate(Detector_t::Type_t::CB | Detector_t::Type_t::PID, 200,2,2,1,0,0,0,{})
        );

  buffer->Candidates.back().Clusters.push_back(
        TCluster(TVector3(25,0,0), 270, -350, Detector_t::Type_t::CB, 0)
        );


  buffer->Candidates.back().Clusters.back().Hits.push_back(
        TClusterHit(110, {
                           {Channel_t::Type_t::Integral, 150}, // MeV
                           {Channel_t::Type_t::Timing, -290}   // ns
                         })
        );

  buffer->Candidates.back().Clusters.back().Hits.push_back(
        TClusterHit(220, {
                           {Channel_t::Type_t::Integral, 120}, // MeV
                           {Channel_t::Type_t::Timing, -280}   // ns
                         })
        );

  buffer->Candidates.back().Clusters.push_back(
        TCluster(TVector3(10,0,0), 5, -270, Detector_t::Type_t::PID, 0)
        );

  buffer->Candidates.back().Clusters.back().Hits.push_back(
        TClusterHit(20, {{Channel_t::Type_t::Integral, 4}})
        );

  //======= Candidate 2 =============
  buffer->Candidates.push_back(
              TCandidate(Detector_t::Type_t::TAPS | Detector_t::Type_t::TAPSVeto, 300,2,3,0.5,0,0,0,{})
        );
  buffer->Candidates.back().Clusters.push_back(
        TCluster(TVector3(25,10,0),300,-200,Detector_t::Type_t::TAPS,0)
        );
  buffer->Candidates.back().Clusters.push_back(
        TCluster(TVector3(10,1,0),5,-250, Detector_t::Type_t::TAPSVeto,0)
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

  TEvent* readback = nullptr;

  tree->SetBranchAddress("event",&readback);

  REQUIRE(tree->GetEntries()==1);

  tree->GetEntry(0);

  cout << *readback << endl;

  /// \todo Compare the events. implement == operator for events...


}
