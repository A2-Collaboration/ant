#include "ProtonTagger.h"
#include "base/std_ext/math.h"

#include "TTree.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;


ProtonTagger::ProtonTagger(const string& name, OptionsPtr opts):
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
        time += c->Time * c->CaloEnergy;
        energy += c->CaloEnergy;
    }
    return time / energy;
}

void ProtonTagger::ProcessEvent(const TEvent& event, manager_t&)
{

    auto taps_cands = event.Reconstructed().Candidates.get_ptr_list(
                          [] (const TCandidate& c) {
        return (c.Detector & Detector_t::Any_t::TAPS_Apparatus) && c.CaloEnergy > 50.0;
    });
    if(taps_cands.size() != 1)
        return;



    TParticleList cb_photons;

    for(const auto& p : event.Reconstructed().Particles.Get(ParticleTypeDatabase::Photon)) {
        if((p->Candidate && p->Candidate->Detector & Detector_t::Any_t::CB_Apparatus) && p->Ek() > 50.0) {
            cb_photons.emplace_back(p);
        }
    }

    if(cb_photons.size() != 2)
        return;



    const auto& gg = *cb_photons.at(0) + *cb_photons.at(1);
    b_ggIM   = gg.M();

    b_cbtime = TimeAverage({cb_photons.at(0)->Candidate,cb_photons.at(1)->Candidate});

    const LorentzVec target(0,0,0,ParticleTypeDatabase::Proton.Mass());

    for(const auto& t : event.Reconstructed().TaggerHits) {

        b_tagTime = t.Time;
        b_tagCh   = t.Channel;

        const ant::LorentzVec mm = t.GetPhotonBeam() + target - gg;

        b_MM = mm.M();

            for(const auto& p : taps_cands) {
                b_Size = p->ClusterSize;
                b_E    = p->CaloEnergy;
                b_veto = p->VetoEnergy;
                b_Time = p->Time;

                b_angle = radian_to_degree(mm.Angle(*p));

                const auto& cluster = p->FindCaloCluster();

                b_Eshort = 0.0;
                b_EshortCentral = 0.0;

                if(cluster) {
                    for(const TClusterHit& clusterhit : cluster->Hits) {
                        for(const TClusterHit::Datum& datum : clusterhit.Data) {

                            if(datum.Type != Channel_t::Type_t::IntegralShort)
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
