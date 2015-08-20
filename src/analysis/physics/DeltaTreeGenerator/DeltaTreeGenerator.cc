#include "DeltaTreeGenerator.h"

#include "plot/root_draw.h"
#include "plot/Histogram.h"
#include "data/Event.h"
#include "data/Particle.h"
#include "base/ParticleType.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;
using namespace std;

DeltaTreeGenerator::DeltaTreeGenerator(const string &name):
    Physics(name)
{
    cout << "DeltaTreeGenerator:" << endl;

    cout << "  Adding Photon-Tree" << endl;
    photonTree    = new TTree("photons","photons");
    reconstructed = new TClonesArray("TLorentzVector",1);
    mctrue        = new TClonesArray("TLorentzVector",1);
    photonTree->Branch("reconstructed",addressof(reconstructed));
    photonTree->Branch("mctrue",addressof(mctrue));


    cout << "  Adding Tagger-Histogram" << endl;
    taggerHits = HistFac.makeTH1D("Tagger Hits","Energy [MeV]","#",BinSettings(200,200,450));
    mcgamma = HistFac.makeTH1D("#gamma - MC true","# #gamma","#",BinSettings(8));
    recgamma = HistFac.makeTH1D("#gamma - reconstructed","# #gamma","#",BinSettings(8));
}

void DeltaTreeGenerator::ProcessEvent(const Event &event)
{
    ParticleList photons;
    ParticleList mc_photons;
    for( const auto& particle : event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon) )
        photons.emplace_back(particle);
    for ( const auto& mc_particle: event.MCTrue().Particles().Get(ParticleTypeDatabase::Photon))
        mc_photons.emplace_back(mc_particle);

    mcgamma->Fill(mc_photons.size());
    recgamma->Fill(photons.size());

    reconstructed->Clear();
    reconstructed->Expand(photons.size());
    for ( size_t i = 0 ; i < photons.size() ; ++i)
    {
        TLorentzVector* lph = (TLorentzVector*)reconstructed->ConstructedAt(i);
        lph->SetPxPyPzE(photons.at(i)->Px(),
                        photons.at(i)->Py(),
                        photons.at(i)->Pz(),
                        photons.at(i)->E());
    }

    mctrue->Clear();
    mctrue->Expand(mc_photons.size());
    for ( size_t i = 0 ; i < mc_photons.size() ; ++i)
    {
        TLorentzVector* lph = (TLorentzVector*)mctrue->ConstructedAt(i);
        lph->SetPxPyPzE(mc_photons.at(i)->Px(),
                        mc_photons.at(i)->Py(),
                        mc_photons.at(i)->Pz(),
                        mc_photons.at(i)->E());
    }
    photonTree->Fill();

    for( const auto& taggerHit: event.MCTrue().TaggerHits())
        taggerHits->Fill(taggerHit->PhotonEnergy());
}

void DeltaTreeGenerator::Finish()
{
}

void DeltaTreeGenerator::ShowResult()
{
    canvas("treegen") << taggerHits << mcgamma << recgamma << endc;
}

AUTO_REGISTER_PHYSICS(DeltaTreeGenerator, "DeltaTreeGenerator")
