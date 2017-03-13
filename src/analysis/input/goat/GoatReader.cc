#include "GoatReader.h"

#include "detail/TreeManager.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/Trigger.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext/memory.h"
#include "base/Detector_t.h"

#include "TTree.h"

#include <string>
#include <iostream>
#include <memory>


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;
using namespace std;

/**
 * @brief map goat apparatus numbers to apparatus_t enum values
 * in case unknown values show up: -> exception and do not sliently ignore
 */
Detector_t::Any_t IntToDetector_t(const int& a) {
    constexpr auto DETECTOR_NaI = 1;
    constexpr auto DETECTOR_PID = 2;
    constexpr auto DETECTOR_MWPC = 4;
    constexpr auto DETECTOR_BaF2 = 8;
    constexpr auto DETECTOR_PbWO4 = 16;
    constexpr auto DETECTOR_Veto = 32;

    auto d = Detector_t::Any_t::None;
    if(a & DETECTOR_NaI) {
        d |= Detector_t::Type_t::CB;
    }
    if(a & DETECTOR_PID) {
        d |= Detector_t::Type_t::PID;
    }
    if(a & DETECTOR_MWPC) {
        d |= Detector_t::Any_t::Tracker;
    }
    if(a & DETECTOR_BaF2) {
        d |= Detector_t::Type_t::TAPS;
    }
    if(a & DETECTOR_PbWO4) {
        d |= Detector_t::Type_t::TAPS;
    }
    if(a & DETECTOR_Veto) {
        d |= Detector_t::Type_t::TAPSVeto;
    }
    return d;
}

void GoatReader::CopyDetectorHits(TEventData& recon)
{
//    auto fill_readhits = [] (TEventData& recon,
//            const DetectorHitInput::Item_t item,
//            Detector_t::Type_t type) {
//        for(int i=0;i<item.nHits;i++) {
//            const unsigned channel = item.Hits[i];
//            const double energy = item.Energy[i];
//            const double time = item.Time[i];
//            if(isfinite(energy)) {
//                recon.DetectorReadHits.emplace_back(
//                            LogicalChannel_t{type, Channel_t::Type_t::Integral, channel},
//                            TDetectorReadHit::Value_t{energy}
//                            );
//            }
//            if(isfinite(time)) {
//                recon.DetectorReadHits.emplace_back(
//                            LogicalChannel_t{type, Channel_t::Type_t::Timing, channel},
//                            TDetectorReadHit::Value_t{time}
//                            );
//            }
//        }
//    };

//    fill_readhits(recon, detectorhits.NaI, Detector_t::Type_t::CB);
//    fill_readhits(recon, detectorhits.PID, Detector_t::Type_t::PID);
//    fill_readhits(recon, detectorhits.BaF2, Detector_t::Type_t::TAPS);
//    fill_readhits(recon, detectorhits.Veto, Detector_t::Type_t::TAPSVeto);
    /// \todo think about MWPC stuff here
}

void GoatReader::CopyTagger(TEventData& recon)
{
//    for( Int_t i=0; i<tagger.GetNTagged(); ++i) {
//        recon.TaggerHits.emplace_back(
//                    tagger.GetTaggedChannel(i),
//                    tagger.GetTaggedEnergy(i),
//                    tagger.GetTaggedTime(i)
//                    );
//    }
}

void GoatReader::CopyTrigger(TEventData& recon)
{
//    TTrigger& ti = recon.Trigger;

//    ti.DAQEventID = eventParameters.EventNumber;
//    ti.CBEnergySum = trigger.GetEnergySum();
//    ti.ClusterMultiplicity = trigger.GetMultiplicity();
//    ti.CBTiming = 0; // Timing is set after event is completely copied

//    for( int err=0; err < trigger.GetNErrors(); ++err) {
//        ti.DAQErrors.emplace_back(
//                    trigger.GetErrorModuleID()[err],
//                    trigger.GetErrorModuleIndex()[err],
//                    trigger.GetErrorCode()[err]
//                    );
//    }
}

/**
 * @brief map the cluster sizes from goat to unsigned ints
 * negative values mean no hit in the calorimeter
 * map those to 0
 */
clustersize_t GoatReader::MapClusterSize(const int& size) {
    return size < 0 ? 0 : size;
}

void GoatReader::CopyTracks(TEventData& recon)
{


//    for(Int_t i=0; i< tracks.GetNTracks(); ++i) {

//        const Detector_t::Any_t det = IntToDetector_t(tracks.GetDetectors(i));

//        // Goat does not provide clusters,
//        // so simulate some with fuzzy logic...
//        TClusterList clusters;
//        /// \todo how does this work with MWPC?

//        if(det & Detector_t::Any_t::Calo) {

//            clusters.emplace_back(
//                        vec3::RThetaPhi(1.0, tracks.GetTheta(i), tracks.GetPhi(i)),
//                        tracks.GetClusterEnergy(i),
//                        tracks.GetTime(i),
//                        det & Detector_t::Type_t::CB ? Detector_t::Type_t::CB : Detector_t::Type_t::TAPS ,
//                        tracks.GetCentralCrystal(i)
//                        );
//            auto& calo_cluster = clusters.back();
//            if(tracks.GetShortEnergy(i)>0)
//                calo_cluster.ShortEnergy = tracks.GetShortEnergy(i);

//            double vetoEnergy = 0.0;
//            if(det & Detector_t::Any_t::Veto) {
//                vetoEnergy =  tracks.GetVetoEnergy(i);
//                clusters.emplace_back(
//                            vec3(std_ext::NaN, std_ext::NaN, std_ext::NaN), // no veto position available
//                            vetoEnergy,
//                            std_ext::NaN, // no veto timing available
//                            det & Detector_t::Type_t::PID ? Detector_t::Type_t::PID : Detector_t::Type_t::TAPSVeto,
//                            tracks.GetCentralVeto(i)
//                            );

//            }

//            recon.Candidates.emplace_back(
//                        det,
//                        tracks.GetClusterEnergy(i),
//                        tracks.GetTheta(i),
//                        tracks.GetPhi(i),
//                        tracks.GetTime(i),
//                        MapClusterSize(tracks.GetClusterSize(i)),
//                        vetoEnergy,
//                        //tracks.GetMWPC0Energy(i)+tracks.GetMWPC1Energy(i),
//                        std_ext::NaN, // MWPC not handled correctly at the moment
//                        TClusterList(clusters.begin(), clusters.end())
//                        );
//        }
//        else if(det & Detector_t::Any_t::Veto) {
//            // veto-only track is just a cluster in Ant
//            // don't know if such tracks actually exist in GoAT/Acqu...
//            const double vetoEnergy =  tracks.GetVetoEnergy(i);
//            clusters.emplace_back(
//                        vec3(std_ext::NaN, std_ext::NaN, std_ext::NaN), // no veto position available
//                        vetoEnergy,
//                        std_ext::NaN, // no veto timing available
//                        det & Detector_t::Type_t::PID ? Detector_t::Type_t::PID : Detector_t::Type_t::TAPSVeto,
//                        tracks.GetCentralVeto(i)
//                        );
//        }

//        // always add clusters...
//        recon.Clusters.insert(recon.Clusters.end(), clusters.begin(), clusters.end());
//    }
}

GoatReader::GoatReader(const std::shared_ptr<const WrapTFileInput>& rootfiles) :
    inputfiles(rootfiles), // remember ROOT files to keep them open as long as WrapTTrees exist
    trigger(ExpConfig::Setup::GetDetector<expconfig::detector::Trigger>())
{
    treeDetectorHitInput.LinkBranches(*inputfiles);
    treeTaggerInput.LinkBranches(*inputfiles);
    treeTriggerInput.LinkBranches(*inputfiles);
    treeTrackInput.LinkBranches(*inputfiles);
}

GoatReader::~GoatReader() {}

bool GoatReader::IsSource() {
//    return trees->GetEntries()>0;
}


bool GoatReader::ReadNextEvent(event_t& event)
{
//    if(current_entry>=trees->GetEntries())
//        return false;



    if(!event.HasReconstructed()) {
        /// \todo think of some better timestamp?
        const TID tid(
                    static_cast<std::uint32_t>(std::time(nullptr)),
                    static_cast<std::uint32_t>(current_entry),
                    std::list<TID::Flags_t>{TID::Flags_t::AdHoc}
                    );
        event.MakeReconstructed(tid);
    }

    auto& recon = event.Reconstructed();

    CopyDetectorHits(recon);
    CopyTrigger(recon);
    CopyTagger(recon);
    CopyTracks(recon);

    // calculate trigger avg timing
    recon.Trigger.CBTiming = trigger->GetCBTiming(recon);

    ++current_entry;
    return true;
}

double GoatReader::PercentDone() const
{
//    return double(current_entry)/double(trees->GetEntries());
}

bool GoatReader::treeDetectorHitInput_t::LinkBranches(const WrapTFileInput& input)
{
    TTree* tree;
    if(!input.GetObject("detectorHits", tree))
        return false;
    NaI.LinkBranches(tree);
    PID.LinkBranches(tree);
    MWPC.LinkBranches(tree);
    BaF2.LinkBranches(tree);
    Veto.LinkBranches(tree);
    return true;
}

bool GoatReader::treeTaggerInput_t::LinkBranches(const WrapTFileInput& input)
{
    if(!input.GetObject("tagger",t.Tree))
        return false;
    t.LinkBranches();
    return true;
}

bool GoatReader::treeTriggerInput_t::LinkBranches(const WrapTFileInput& input)
{
    if(!input.GetObject("eventParameters",tEventParams.Tree))
        return false;
    if(!input.GetObject("trigger",t.Tree))
        return false;
    tEventParams.LinkBranches();
    t.LinkBranches();
    return true;
}

bool GoatReader::treeTrackInput_t::LinkBranches(const WrapTFileInput& input)
{
    if(!input.GetObject("tracks", t.Tree))
        return false;
    t.LinkBranches();
    return true;
}
