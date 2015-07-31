#include "TAPS_Energy.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDataRecord.h"

#include <list>


using namespace std;
using namespace ant;
using namespace ant::calibration;

TAPS_Energy::TAPS_Energy(std::shared_ptr<CalibrationDataManager> calmgr,
                         Calibration::Converter::ptr_t converter,
                         double defaultPedestal,
                         double defaultGain,
                         double defaultThreshold,
                         double defaultRelativeGain):
    Energy(Detector_t::Type_t::TAPS,
           calmgr,
           converter,
           defaultPedestal,
           defaultGain,
           defaultThreshold,
           defaultRelativeGain)
{

}

TAPS_Energy::ThePhysics::ThePhysics(const string& name):
    Physics(name)
{
    const BinSettings cb_channels(438);
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH2D("2 neutral IM (TAPS,CB)", "IM [MeV]", "#", energybins, cb_channels, "ggIM");
}

void TAPS_Energy::ThePhysics::ProcessEvent(const Event& event)
{
    const auto& cands = event.Reconstructed().Candidates();

    const auto CBTAPS = detector_t::CB | detector_t::TAPS;

    for( auto comb = makeCombination(cands,2); !comb.Done(); ++comb ) {
        const CandidatePtr& p1 = comb.at(0);
        const CandidatePtr& p2 = comb.at(1);

        if(p1->VetoEnergy()==0 && p2->VetoEnergy()==0) {

            //require exactly 1 CB and 1 TAPS
            const detector_t dets = (p1->Detector() & CBTAPS) ^ (p2->Detector() & CBTAPS);

            if(dets & CBTAPS) {
                const Particle a(ParticleTypeDatabase::Photon,comb.at(0));
                const Particle b(ParticleTypeDatabase::Photon,comb.at(1));
                const TLorentzVector gg = a + b;

                // Find the one that was in TAPS
                auto cl = p1->Detector() & detector_t::TAPS ? p1->FindCaloCluster() : p2->FindCaloCluster();
                if(cl)
                    ggIM->Fill(gg.M(),cl->CentralElement);
            }
        }
    }

}

void TAPS_Energy::ThePhysics::Finish()
{
}

void TAPS_Energy::ThePhysics::ShowResult()
{
    canvas("TAPS_Energy") << drawoption("colz") << ggIM << endc;
}

unique_ptr<Physics> TAPS_Energy::GetPhysicsModule()
{
    return std_ext::make_unique<ThePhysics>(GetName());
}
