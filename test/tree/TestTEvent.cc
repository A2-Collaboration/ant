#include "catch.hpp"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "base/tmpfile_t.h"
#include "base/std_ext/memory.h"

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
  auto event = new TEvent();

  tree->Branch(branchname.c_str(), event);

  event->MakeReconstructed(TID(10));
  auto& eventdata = event->Reconstructed();

  eventdata.DetectorReadHits.emplace_back();
  eventdata.DetectorReadHits.emplace_back();
  eventdata.DetectorReadHits.emplace_back();

  auto& clusters = eventdata.Clusters;

  clusters.emplace_back(vec3(1,2,3),
                        100, 0.5,
                        Detector_t::Type_t::PID,
                        127, // central element
                        vector<TClusterHit>{TClusterHit()}
                        );
  clusters.emplace_back(vec3(4,5,6),
                        100, 0.5,
                        Detector_t::Type_t::CB,
                        127, // central element
                        vector<TClusterHit>{TClusterHit(), TClusterHit()}
                        );
  clusters.emplace_back(vec3(7,8,9),
                        100, 0.5,
                        Detector_t::Type_t::TAPS,
                        127, // central element
                        vector<TClusterHit>{TClusterHit(), TClusterHit(), TClusterHit()}
                        );

  auto cluster0 = std::next(clusters.begin(), 0);
  auto cluster1 = std::next(clusters.begin(), 1);
  auto cluster2 = std::next(clusters.begin(), 2);

  auto& candidates = eventdata.Candidates;

  candidates.emplace_back(
              Detector_t::Any_t::CB_Apparatus,
              200,
              0.0, 0.0, 0.0, // theta/phi/time
              2, // cluster size
              2.0, 0.0, // veto/tracker
              TClusterList{cluster1, cluster0}
              );
  candidates.emplace_back(
              Detector_t::Any_t::TAPS_Apparatus,
              100,
              1.0, 2.0, 3.0, // theta/phi/time
              8, // cluster size
              2.0, 0.0, // veto/tracker
              TClusterList{cluster2}
              );

  auto candidate0 = std::next(candidates.begin(), 0);

  auto particle0 = make_shared<TParticle>(ParticleTypeDatabase::Photon, candidate0.get_ptr());
  auto particle1 = make_shared<TParticle>(ParticleTypeDatabase::Photon, LorentzVec(7,8,9,10));
  auto particle2 = make_shared<TParticle>(ParticleTypeDatabase::Pi0, LorentzVec(3,4,5,6));


  eventdata.Particles.Add(particle0);
  eventdata.Particles.Add(particle1);

  eventdata.ParticleTree = Tree<TParticlePtr>::MakeNode(particle2);
  eventdata.ParticleTree->CreateDaughter(particle1);
  eventdata.ParticleTree->CreateDaughter(particle0);

  cout << event << endl;
  cout << *event << endl;

  REQUIRE(event->Reconstructed().Candidates.size()==2);

  tree->Fill();

  event->MakeReconstructed(TID());
  event->MakeMCTrue(TID());

  event->Reconstructed().Particles.Add(particle0);
  event->Reconstructed().Particles.Add(particle0);

  event->MCTrue().Particles.Add(particle0);
  event->MCTrue().Particles.Add(particle1);

  tree->Fill();

  f.Write();
  f.Close();

  tree = nullptr;

  TFile f2(tmpfile.filename.c_str(),"READ");
  REQUIRE(f2.IsOpen());

  f2.GetObject(treename.c_str(), tree);

  REQUIRE(tree!=nullptr);

  TEvent* readback_event = nullptr;

  tree->SetBranchAddress(branchname.c_str(), &readback_event);

  REQUIRE(tree->GetEntries() == 2);

  tree->GetEntry(0);

  cout << readback_event << endl;
  cout << *readback_event << endl;


  REQUIRE(readback_event->HasReconstructed());
  const auto& readback = readback_event->Reconstructed();


  REQUIRE(readback.ID == TID(10));

  REQUIRE(readback.DetectorReadHits.size() == 3);

  REQUIRE(readback.Clusters.size() == 3);
  REQUIRE(readback.Clusters.at(0).Position == vec3(1,2,3));
  REQUIRE(readback.Clusters.at(2).Position == vec3(7,8,9));
  REQUIRE(readback.Clusters.at(0).Hits.size() == 1);
  REQUIRE(readback.Clusters.at(2).Hits.size() == 3);

  REQUIRE(readback.Candidates.size() == 2);

  REQUIRE(readback.Clusters.get_ptr_at(0) == readback.Candidates.at(0).Clusters.get_ptr_at(1));

  REQUIRE(readback.Particles.GetAll().size() == 2);
  REQUIRE(readback.Particles.GetAll().at(0)->Type() == ParticleTypeDatabase::Photon);
  REQUIRE(readback.Particles.GetAll().at(0)->Candidate == readback.Candidates.get_ptr_at(0));
  REQUIRE(readback.Particles.GetAll().at(1)->E == 10);
  REQUIRE(readback.Particles.Get(ParticleTypeDatabase::Photon).size() == 2);
  REQUIRE(readback.Particles.Get(ParticleTypeDatabase::Proton).size() == 0);

  REQUIRE(readback.ParticleTree != nullptr);
  // check if that particle was properly re-created
  REQUIRE(readback.ParticleTree->Get() != particle2);
  REQUIRE(readback.ParticleTree->Get()->Type() == particle2->Type());
  REQUIRE(readback.ParticleTree->Get()->Type() == ParticleTypeDatabase::Pi0);
  REQUIRE(readback.ParticleTree->Daughters().size() == 2);
  REQUIRE(readback.ParticleTree->Daughters().front()->Get() == readback.Particles.GetAll().at(1));
  REQUIRE(readback.ParticleTree->Daughters().back()->Get() == readback.Particles.GetAll().at(0));


  // check some list capabilities

  TCandidatePtrList list;
  for(auto it_cand : readback.Candidates.get_iter()) {
      REQUIRE(it_cand.get_ptr());
      REQUIRE_FALSE(it_cand->Clusters.empty());
      list.emplace_back(it_cand);
  }
  REQUIRE(list.size() == 2);

  auto all_cands = readback.Candidates.get_ptr_list();
  REQUIRE(all_cands.size() == 2);

  auto taps_cands = readback.Candidates.get_ptr_list(
                        [] (const TCandidate& c) { return c.Detector & Detector_t::Type_t::TAPS; } );
  REQUIRE(taps_cands.size() == 1);

}
