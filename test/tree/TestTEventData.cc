#include "catch.hpp"

#include "tree/TEventData.h"

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

  const std::string treename = "t";
  const std::string branchname = "b";

  TFile f(tmpfile.filename.c_str(),"RECREATE");

  TTree* tree = new TTree(treename.c_str(),"");
  TEventData* eventdata = new TEventData();

  tree->Branch(branchname.c_str(), eventdata);

  eventdata->ID = TID(10);

  eventdata->DetectorHits.emplace_back();
  eventdata->DetectorHits.emplace_back();
  eventdata->DetectorHits.emplace_back();

  auto cluster1 = make_shared<TCluster>(TVector3(1,2,3),
                                       100, 0.5,
                                       Detector_t::Type_t::PID,
                                       127, // central element
                                       vector<TClusterHit>{TClusterHit()}
                                       );
  auto cluster2 = make_shared<TCluster>(TVector3(4,5,6),
                                       100, 0.5,
                                       Detector_t::Type_t::CB,
                                       127, // central element
                                       vector<TClusterHit>{TClusterHit(), TClusterHit()}
                                       );
  auto cluster3 = make_shared<TCluster>(TVector3(7,8,9),
                                       100, 0.5,
                                       Detector_t::Type_t::TAPS,
                                       127, // central element
                                       vector<TClusterHit>{TClusterHit(), TClusterHit(), TClusterHit()}
                                       );
  eventdata->Clusters.emplace_back(cluster1);
  eventdata->Clusters.emplace_back(cluster2);
  eventdata->Clusters.emplace_back(cluster3);



  eventdata->Candidates.emplace_back(make_shared<TCandidate>(
                                         Detector_t::Any_t::CB_Apparatus,
                                         200,
                                         0.0, 0.0, 0.0, // theta/phi/time
                                         2, // cluster size
                                         2.0, 0.0, // veto/tracker
                                         vector<TClusterPtr>{cluster2, cluster1}
                                         ));




  tree->Fill();

  cout << eventdata << endl;
  cout << *eventdata << endl;

  f.Write();
  f.Close();

  tree = nullptr;

  TFile f2(tmpfile.filename.c_str(),"READ");
  REQUIRE(f2.IsOpen());

  f2.GetObject(treename.c_str(), tree);

  REQUIRE(tree!=nullptr);

  TEventData* readback = nullptr;

  tree->SetBranchAddress(branchname.c_str(), &readback);

  REQUIRE(tree->GetEntries()==1);

  tree->GetEntry(0);

  cout << readback << endl;
  cout << *readback << endl;

  REQUIRE(readback->ID == TID(10));

  REQUIRE(readback->DetectorHits.size() == 3);

  REQUIRE(readback->Clusters.size() == 3);
  REQUIRE(readback->Clusters.at(0)->Position == TVector3(1,2,3));
  REQUIRE(readback->Clusters.at(2)->Position == TVector3(7,8,9));
  REQUIRE(readback->Clusters.at(0)->Hits.size() == 1);
  REQUIRE(readback->Clusters.at(2)->Hits.size() == 3);

  REQUIRE(readback->Candidates.size() == 1);

  REQUIRE(readback->Clusters.at(0) == readback->Candidates.at(0)->Clusters.at(1));




}
