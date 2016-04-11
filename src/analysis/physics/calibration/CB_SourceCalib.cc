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
    const BinSettings ADCBins(140);

    HitsADC = HistFac.makeTH2D("Detector Hits for ADC-Channel", "ADC-Channel", "#",
                            ADCBins, cb_channels, "HitsADC");
    h_cbdisplay = HistFac.make<TH2CB>("h_cbdisplay","Number of entries");
    KristallHits= HistFac.makeTH1D("Hits in one Crystall", "ADC-Channel", "Hits", ADCBins, "KristallHits");

}
void CB_SourceCalib::ProcessEvent(const TEvent& event, manager_t&)

{

    for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
        // ignore readhits which are not from CB
        if(readhit.DetectorType != Detector_t::Type_t::CB)
            continue;

        //cout << readhit << endl;
        for( int Channel=1; Channel<=720; ++Channel)
            if (readhit.Channel==Channel && readhit.ChannelType == Channel_t::Type_t::Integral)
            {
                const auto& values = adc_converter.Convert(readhit.RawData);
                KristallHits->Fill(values.at(0));
                ;
            }

//       if (readhit.Channel==483 && readhit.ChannelType == Channel_t::Type_t::Integral )
//      {
//           const auto& values = adc_converter.Convert(readhit.RawData);
//           KristallHits->Fill(values.at(0));
//       }
        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            const auto& values = adc_converter.Convert(readhit.RawData);
            cout << values.at(0) << endl;
            HitsADC->Fill(values.at(0),readhit.Channel);
        }
    }




}


void CB_SourceCalib::ShowResult()
{
   h_cbdisplay->SetElements(*HitsADC->ProjectionY());
   canvas(GetName()) << drawoption("colz") << HitsADC
                      << KristallHits
                      << h_cbdisplay
                      << endc;
}

AUTO_REGISTER_PHYSICS(CB_SourceCalib)
