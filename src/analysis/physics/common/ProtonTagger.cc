#include "ProtonTagger.h"
#include "base/std_ext/math.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;


ProtonTagger::ProtonTagger(const string& name, PhysOptPtr opts):
    Physics(name, opts)
{
    const BinSettings binE(1000);
    const BinSettings binVeto(20,0,20);
    const BinSettings binT(1000,-20,20);
    const BinSettings binMM(1000, 500, 1500);
    const BinSettings binAngle(180,0,90.0);

    tof  = HistFac.makeTH2D("ToF", "Time [ns]", "Energy [MeV]", binT, binE, "tof");
    dEE  = HistFac.makeTH2D("dEE", "Energy [MeV]", "Veto Energy [MeV]", binE, binVeto, "dEE");
    cls  = HistFac.makeTH2D("Cluster Size", "Energy [MeV]", "ClusterSize", binE, BinSettings(20), "clusterSize");


    ggIM = HistFac.makeTH1D("2#gamma IM (CB)", "2#gamma IM [MeV]", "", binE, "ggIM");
    MM_after_cut = HistFac.makeTH1D("MM of 2#gamma (CB)", "MM [MeV]", "", binMM, "MM");
    angle = HistFac.makeTH1D("Angle", "[#circ]", "", binAngle, "angle");
}

void ProtonTagger::ProcessEvent(const data::Event& event)
{

    data::CandidateList taps_hits;

    for(const auto& p : event.Reconstructed().Candidates()) {
        if(p->Detector() & Detector_t::Any_t::TAPS) {
            taps_hits.emplace_back(p);
        }
    }
    if(taps_hits.size() != 1)
        return;

    data::ParticleList cb_photons;

    for(const auto& p : event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon)) {
        if(p->Candidate() && p->Candidate()->Detector() & Detector_t::Any_t::CB) {
            cb_photons.emplace_back(p);
        }
    }

    if(cb_photons.size() != 2)
        return;



    const TLorentzVector gg = *cb_photons.at(0) + *cb_photons.at(1);
    ggIM->Fill(gg.M());

    if(!pi0_cut.Contains(gg.M()))
        return;

    const TLorentzVector target(0,0,0,ParticleTypeDatabase::Proton.Mass());

    for(const auto& t : event.Reconstructed().TaggerHits()) {
        const TLorentzVector mm = t->PhotonBeam() + target - gg;
        MM_after_cut->Fill(mm.M());
        if(mm_cut.Contains(mm.M())) {
            for(const auto& p : taps_hits) {

                const auto a = radian_to_degree(mm.Angle(*p));

                angle->Fill(a);

                if(a < 10.0) {
                    tof->Fill(p->Time(), p->ClusterEnergy());
                    dEE->Fill(p->ClusterEnergy(), p->VetoEnergy());
                    cls->Fill(p->ClusterEnergy(), p->ClusterSize());
                }
            }
        }
    }

}

void ProtonTagger::ShowResult()
{
    canvas("Proton Tagger") << drawoption("colz") << tof << dEE << cls << ggIM << MM_after_cut << angle << endc;

}

AUTO_REGISTER_PHYSICS(ProtonTagger)
