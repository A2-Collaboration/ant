#include "CB_SourceCalib.h"

#include "utils/combinatorics.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"


using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

CB_SourceCalib::CB_SourceCalib(const string & name, OptionsPtr opts) :
    Physics(name, opts)
{
    auto detector = ExpConfig::Setup::GetDetector(Detector_t::Type_t::CB);

    const BinSettings cb_channels(detector->GetNChannels());
    const BinSettings energybins(1000);

    ggIM = HistFac.makeTH2D("2 neutral IM (CB,CB)", "IM [MeV]", "#",
                            energybins, cb_channels, "ggIM");
    h_cbdisplay = HistFac.make<TH2CB>("h_cbdisplay","Number of entries");


}
void CB_SourceCalib::ProcessEvent(const TEvent& event, manager_t&)

{

    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        // ignore readhits which are not from CB
        if(readhit.DetectorType != Detector_t::Type_t::CB)
            continue;

        cout << readhit << endl;

        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            const auto& values = adc_converter.Convert(readhit.RawData);
            cout << values.at(0) << endl;
            ggIM->Fill(values.at(0),readhit.Channel);
        }
    }
//        if(readhit.DetectorType != Detector->Type)
//            continue;
//        if(readhit.Values.empty())
//            continue;
//        auto& hit = hits[readhit.Channel];
//        if(readhit.ChannelType == Channel_t::Type_t::Integral)
//            hit.Energy = readhit.Values.front();
//        if(readhit.ChannelType == Channel_t::Type_t::Timing)
//            hit.Time = readhit.Values.front();

    //const auto& cands = event.Reconstructed().Candidates;


//    for( auto comb = analysis::utils::makeCombination(cands.get_ptr_list(),2); !comb.Done(); ++comb ) {
//        const TCandidatePtr& p1 = comb.at(0);
//        const TCandidatePtr& p2 = comb.at(1);

//        if(p1->VetoEnergy==0   && p2->VetoEnergy==0
//           && (p1->Detector & Detector_t::Type_t::CB)
//           && (p2->Detector & Detector_t::Type_t::CB)) {
//            const TParticle a(ParticleTypeDatabase::Photon,comb.at(0));
//            const TParticle b(ParticleTypeDatabase::Photon,comb.at(1));
//            const auto& gg = a+b;

//            auto cl1 = p1->FindCaloCluster();
//            if(cl1)
//                ggIM->Fill(gg.P(),cl1->CentralElement);

//            auto cl2 = p2->FindCaloCluster();
//            if(cl2)
//                ggIM->Fill(gg.P(),cl2->CentralElement);
//        }
//    }
}


void CB_SourceCalib::ShowResult()
{
   h_cbdisplay->SetElements(*ggIM->ProjectionY());
   canvas(GetName()) << drawoption("colz") << ggIM
                      << h_cbdisplay
                      << endc;
}

AUTO_REGISTER_PHYSICS(CB_SourceCalib)
