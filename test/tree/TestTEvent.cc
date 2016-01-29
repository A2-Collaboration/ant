#include "catch.hpp"

#include "tree/TEvent.h"

#include "base/tmpfile_t.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("TEvent: Write/Read TTree", "[tree]") {
  dotest();
}

void dotest() {
  tmpfile_t tmpfile;

  const std::string treename = "t";
  const std::string branchname = "b";

  TFile f(tmpfile.filename.c_str(),"RECREATE");

  TTree* tree = new TTree(treename.c_str(),"");
  TEvent* event = new TEvent();

  tree->Branch(branchname.c_str(), event);

  event->Reconstructed = make_shared<TEventData>();
  auto& eventdata = event->Reconstructed;

  eventdata->ID = TID(10);

  eventdata->DetectorHits.emplace_back();
  eventdata->DetectorHits.emplace_back();
  eventdata->DetectorHits.emplace_back();

  auto cluster0 = make_shared<TCluster>(TVector3(1,2,3),
                                       100, 0.5,
                                       Detector_t::Type_t::PID,
                                       127, // central element
                                       vector<TClusterHit>{TClusterHit()}
                                       );
  auto cluster1 = make_shared<TCluster>(TVector3(4,5,6),
                                       100, 0.5,
                                       Detector_t::Type_t::CB,
                                       127, // central element
                                       vector<TClusterHit>{TClusterHit(), TClusterHit()}
                                       );
  auto cluster2 = make_shared<TCluster>(TVector3(7,8,9),
                                       100, 0.5,
                                       Detector_t::Type_t::TAPS,
                                       127, // central element
                                       vector<TClusterHit>{TClusterHit(), TClusterHit(), TClusterHit()}
                                       );
  eventdata->Clusters.emplace_back(cluster0);
  eventdata->Clusters.emplace_back(cluster1);
  eventdata->Clusters.emplace_back(cluster2);


  auto candidate0 = make_shared<TCandidate>(
                        Detector_t::Any_t::CB_Apparatus,
                        200,
                        0.0, 0.0, 0.0, // theta/phi/time
                        2, // cluster size
                        2.0, 0.0, // veto/tracker
                        vector<TClusterPtr>{cluster1, cluster0}
                        );

  eventdata->Candidates.emplace_back(candidate0);

  auto particle0 = make_shared<TParticle>(ParticleTypeDatabase::Photon, candidate0);
  auto particle1 = make_shared<TParticle>(ParticleTypeDatabase::Photon, TLorentzVector(7,8,9,10));
  auto particle2 = make_shared<TParticle>(ParticleTypeDatabase::Pi0, TLorentzVector(3,4,5,6));


  eventdata->Particles.emplace_back(particle0);
  eventdata->Particles.emplace_back(particle1);

  eventdata->ParticleTree = Tree<TParticlePtr>::MakeNode(particle2);
  eventdata->ParticleTree->CreateDaughter(particle1);
  eventdata->ParticleTree->CreateDaughter(particle0);


  tree->Fill();

  cout << event << endl;
  cout << *event << endl;

  f.Write();
  f.Close();

  tree = nullptr;

  TFile f2(tmpfile.filename.c_str(),"READ");
  REQUIRE(f2.IsOpen());

  f2.GetObject(treename.c_str(), tree);

  REQUIRE(tree!=nullptr);

  TEvent* readback_event = nullptr;

  tree->SetBranchAddress(branchname.c_str(), &readback_event);

  REQUIRE(tree->GetEntries()==1);

  tree->GetEntry(0);

  cout << readback_event << endl;
  cout << *readback_event << endl;


  auto& readback = readback_event->Reconstructed;
  REQUIRE(readback != nullptr);


  REQUIRE(readback->ID == TID(10));

  REQUIRE(readback->DetectorHits.size() == 3);

  REQUIRE(readback->Clusters.size() == 3);
  REQUIRE(readback->Clusters.at(0)->Position == TVector3(1,2,3));
  REQUIRE(readback->Clusters.at(2)->Position == TVector3(7,8,9));
  REQUIRE(readback->Clusters.at(0)->Hits.size() == 1);
  REQUIRE(readback->Clusters.at(2)->Hits.size() == 3);

  REQUIRE(readback->Candidates.size() == 1);

  REQUIRE(readback->Clusters.at(0) == readback->Candidates.at(0)->Clusters.at(1));

  REQUIRE(readback->Particles.size() == 2);
  REQUIRE(readback->Particles.at(0)->Type() == ParticleTypeDatabase::Photon);
  REQUIRE(readback->Particles.at(0)->Candidate == readback->Candidates.at(0));

  REQUIRE(readback->ParticleTree != nullptr);
  // check if that particle was properly re-created
  REQUIRE(readback->ParticleTree->Get() != particle2);
  REQUIRE(readback->ParticleTree->Get()->Type() == particle2->Type());
  REQUIRE(readback->ParticleTree->Get()->Type() == ParticleTypeDatabase::Pi0);
  REQUIRE(readback->ParticleTree->Daughters().size() == 2);
  REQUIRE(readback->ParticleTree->Daughters().front()->Get() == readback->Particles.at(1));
  REQUIRE(readback->ParticleTree->Daughters().back()->Get() == readback->Particles.at(0));




}
