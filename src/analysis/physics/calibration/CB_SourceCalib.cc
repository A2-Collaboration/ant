#include "CB_SourceCalib.h"

#include "utils/combinatorics.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/Trigger.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::physics;

CB_SourceCalib::CB_SourceCalib(const string & name, OptionsPtr opts) :
    Physics(name, opts),
    tdc_converter(expconfig::detector::Trigger::Reference_CATCH_CBCrate)
{
    cb_detector = ExpConfig::Setup::GetDetector<expconfig::detector::CB>();

    const BinSettings cb_channels(cb_detector->GetNChannels());
    const BinSettings ADCBins(140);

    HitsADC = HistFac.makeTH2D("Detector Hits for ADC-Channel", "ADC-Channel", "#",
                            ADCBins, cb_channels, "HitsADC");
    h_cbdisplay = HistFac.make<TH2CB>("h_cbdisplay","Number of entries");
    KristallHits= HistFac.makeTH1D("Hits in one Crystall", "ADC-Channel", "Hits", ADCBins, "KristallHits");
    TimeHits=HistFac.makeTH2D("Time Hits for ADC-Channel", "ADC-Channel", "#", 500, cb_channels, "Time in ns");
    HitsADC_Cluster = HistFac.makeTH2D("Detector Hits for ADC-Channel", "ADC-Channel", "#",
                            ADCBins, cb_channels, "HitsADC_Cluster");
    VerworfeneHits = HistFac.makeTH2D(" Verworfene Einträge", "ADC-Channel", "#", ADCBins, cb_channels, "VerworfeneHits");

}
void CB_SourceCalib::ProcessEvent(const TEvent& event, manager_t&)

{

    // build sorted readhits with refhit only
//    {
//        vector<TDetectorReadHit> readhits;
//        ReconstructHook::Base::readhits_t sorted_readhits;
//        for(const TDetectorReadHit& readhit : event.Reconstructed().DetectorReadHits) {
//            auto& cb_reftiming = expconfig::detector::Trigger::Reference_CATCH_CBCrate.LogicalChannel;
//            if(readhit.DetectorType == cb_reftiming.DetectorType) {
//                readhits.push_back(readhit);
//                sorted_readhits.add_item(readhit.DetectorType, readhits.back());
//            }
//        }
//        tdc_converter.ApplyTo(sorted_readhits);

//    }

    auto& readhits = event.Reconstructed().DetectorReadHits;
    auto& clusters = event.Reconstructed().Clusters;
    // auto cb_element = cb_detector->GetClusterElement(5);
    // cout << "Nachbarn von Cluster Element 5" << cb_element->Neighbours << endl;
    for (const TCluster& cluster : clusters)
    {
        if(cluster.DetectorType != Detector_t::Type_t::CB )
            continue;

        //cout << cluster.Time << endl;
        //cout << cluster.CentralElement <<endl;
//                for(const TClusterHit& hit : cluster.Hits)
//                {

//                    cout<<hit.Time<< "Hit time" << endl;
//                    HitsADC_Cluster->Fill(hit.Energy,hit.Channel);

//                }
        //cout << cluster.Hits << "ClusterHits" << endl;

        if (cluster.Hits.size()==1)
        {
            //cout << cluster.Hits << "  CLusterTime" << endl;
            HitsADC_Cluster->Fill(cluster.Energy,cluster.CentralElement);
            TimeHits->Fill(cluster.Time, cluster.CentralElement);
        }
        else
        {
            VerworfeneHits->Fill(cluster.Energy, cluster.CentralElement);
        }
    }


    for(const TDetectorReadHit& readhit : readhits) {
        // ignore readhits which are not from CB
        if(readhit.DetectorType != Detector_t::Type_t::CB)
            continue;


//        if(readhit.ChannelType == Channel_t::Type_t::Timing) {
//            const auto& timings = tdc_converter.Convert(readhit.RawData);
//            cout << timings << endl;
//        }
            // plot the sum of all hits of the cb aginst the ADC-channel
//          if (readhit.ChannelType == Channel_t::Type_t::Integral)
//            {

//                const auto& values = adc_converter.Convert(readhit.RawData);
//                KristallHits->Fill(values.at(0));


//            }

//       if (readhit.Channel==483 && readhit.ChannelType == Channel_t::Type_t::Integral )
//      {
//           const auto& values = adc_converter.Convert(readhit.RawData);
//           KristallHits->Fill(values.at(0));
//       }

        // plot the Detector elements against the ADC Channel
        if(readhit.ChannelType == Channel_t::Type_t::Integral) {
            const auto& values = adc_converter.Convert(readhit.RawData);
            //cout << values.at(0) << endl;
            HitsADC->Fill(values.at(0),readhit.Channel);
        }
    }




}


void CB_SourceCalib::ShowResult()
{
   h_cbdisplay->SetElements(*HitsADC->ProjectionY());
   canvas(GetName()) << drawoption("colz") << HitsADC
                      << KristallHits
                      << VerworfeneHits
                      << TimeHits
                      << HitsADC_Cluster
                      << endc;
}

AUTO_REGISTER_PHYSICS(CB_SourceCalib)
