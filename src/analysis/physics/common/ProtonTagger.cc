#include "ProtonTagger.h"
#include "base/std_ext/math.h"

#include "TTree.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace ant::analysis::physics;
using namespace std;


ProtonTagger::ProtonTagger(const string& name, PhysOptPtr opts):
    Physics(name, opts)
{

    tree = HistFac.makeTTree("tree");

    tree->Branch("tagTime", &b_tagTime);
    tree->Branch("tagCh",   &b_tagCh);
    tree->Branch("ggIM",    &b_ggIM);
    tree->Branch("MM",      &b_MM);
    tree->Branch("angle",   &b_angle);
    tree->Branch("time",    &b_Time);
    tree->Branch("E",       &b_E);
    tree->Branch("veto",    &b_veto);
    tree->Branch("size",    &b_Size);
    tree->Branch("cbtime",  &b_cbtime);
    tree->Branch("Eshort",  &b_Eshort);
    tree->Branch("EshortCentral",  &b_EshortCentral);
}

template <typename T>
double TimeAverage(const std::initializer_list<const T> cands) {
    double time   = 0.0;
    double energy = 0.0;
    for(const auto& c : cands) {
        time += c->Time * c->ClusterEnergy;
        energy += c->ClusterEnergy;
    }
    return time / energy;
}

void ProtonTagger::ProcessEvent(const data::Event& event)
{

    data::CandidateList taps_hits;

    for(const auto& p : event.Reconstructed.Candidates) {
        if((p->Detector & Detector_t::Any_t::TAPS) && p->ClusterEnergy > 50.0) {
            taps_hits.emplace_back(p);
        }
    }
    if(taps_hits.size() != 1)
        return;



    data::ParticleList cb_photons;

    for(const auto& p : event.Reconstructed.Particles.Get(ParticleTypeDatabase::Photon)) {
        if((p->Candidate && p->Candidate->Detector & Detector_t::Any_t::CB) && p->Energy() > 50.0) {
            cb_photons.emplace_back(p);
        }
    }

    if(cb_photons.size() != 2)
        return;



    const TLorentzVector gg = *cb_photons.at(0) + *cb_photons.at(1);
    b_ggIM   = gg.M();

    b_cbtime = TimeAverage({cb_photons.at(0)->Candidate,cb_photons.at(1)->Candidate});

    const TLorentzVector target(0,0,0,ParticleTypeDatabase::Proton.Mass());

    for(const auto& t : event.Reconstructed.TaggerHits) {

        b_tagTime = t.Time;
        b_tagCh   = t.Channel;

        const TLorentzVector mm = t.GetPhotonBeam() + target - gg;

        b_MM = mm.M();

            for(const auto& p : taps_hits) {
                b_Size = p->ClusterSize;
                b_E    = p->ClusterEnergy;
                b_veto = p->VetoEnergy;
                b_Time = p->Time;

                b_angle = radian_to_degree(mm.Angle(*p));

                const auto& cluster = p->FindCaloCluster();

                b_Eshort = 0.0;
                b_EshortCentral = 0.0;

                if(cluster) {
                    for(const Cluster::Hit& clusterhit : cluster->Hits) {
                        for(const Cluster::Hit::Datum& datum : clusterhit.Data) {

                            if(datum.Type != Channel_t::Type_t::PedestalShort)
                                continue;

                            b_Eshort += datum.Value;
                            if(clusterhit.Channel == cluster->CentralElement) {
                                b_EshortCentral = datum.Value;
                            }
                        }
                    }
                }

                tree->Fill();

            }
    }

}

void ProtonTagger::ShowResult()
{

}

AUTO_REGISTER_PHYSICS(ProtonTagger)
