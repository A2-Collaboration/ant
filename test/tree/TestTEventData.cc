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
  auto cluster1 = make_shared<TCluster>(TVector3(1,2,3),
                                       100, 0.5,
                                       Detector_t::Type_t::CB,
                                       127, // central element
                                       vector<TClusterHit>{TClusterHit()}
                                       );
  auto cluster2 = make_shared<TCluster>(TVector3(1,2,3),
                                       100, 0.5,
                                       Detector_t::Type_t::CB,
                                       127, // central element
                                       vector<TClusterHit>{TClusterHit()}
                                       );
  eventdata->Clusters.emplace_back(cluster1);
  eventdata->Clusters.emplace_back(cluster1);
  eventdata->Clusters.emplace_back(cluster1);
  eventdata->Clusters.emplace_back(cluster2);



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
  REQUIRE(readback->Clusters.size() == 4);
  REQUIRE(readback->Clusters.back()->Position == TVector3(1,2,3));
  REQUIRE(readback->Clusters.back()->Hits.size() == 1);




}
