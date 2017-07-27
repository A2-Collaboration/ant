#include "GoatReader.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "expconfig/ExpConfig.h"
#include "expconfig/detectors/Trigger.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/Detector_t.h"

#include "TTree.h"

#include <string>


using namespace ant;
using namespace ant::analysis::input;
using namespace std;

// for GoatReader::trees std::set
bool operator<(const reference_wrapper<TTree>& t1, const reference_wrapper<TTree>& t2) {
    return addressof(t1.get()) < addressof(t2.get());
}

GoatReader::GoatReader(const std::shared_ptr<const WrapTFileInput>& rootfiles) :
    current_entry(0),
    max_entries(0),
    init(true)
{
    // let those components to the work, collect trees
    init &= treeDetectorHitInput.LinkBranches(*rootfiles, trees);
    init &= treeTaggerInput.LinkBranches(*rootfiles, trees);
    init &= treeTriggerInput.LinkBranches(*rootfiles, trees);
    init &= treeTrackInput.LinkBranches(*rootfiles, trees);

    // return silently if we haven't initiliazed
    if(!init)
        return;

    if(trees.empty())
        throw Exception("No trees were linked");

    // look for trigger (might fail) only after we're sure we match to the file
    inputfiles = rootfiles; // remember ROOT files to keep them open as long as WrapTTrees exist
    trigger = ExpConfig::Setup::GetDetector<expconfig::detector::Trigger>();


    // check tree entries pairwise, and gather unique_TTrees for later usage
    {
        auto it = trees.begin();
        auto it_next = std::next(it);
        for( ; it_next != trees.end() ; ++it, ++it_next) {
            auto& t1 = it->get();
            auto& t2 = it_next->get();
            if(t1.GetEntries() != t2.GetEntries())
                throw Exception("Tree "+string(t1.GetName())+ " and tree "+
                                string(t2.GetName())+" do not have equal entries");
        }
    }

    // as all trees have same number of entries, max_entries is given by front element
    max_entries = trees.begin()->get().GetEntries();

    LOG(INFO) << "Successfully opened GoAT file with " << max_entries << " entries";
}

GoatReader::~GoatReader() {}

reader_flags_t GoatReader::GetFlags() const {
    if(init)
        return reader_flag_t::IsSource;
    else
        return {};
}


bool GoatReader::ReadNextEvent(event_t& event)
{
    if(current_entry>=max_entries)
        return false;

    for(auto& t : trees) {
        t.get().GetEntry(current_entry);
    }

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

    treeDetectorHitInput.Copy(recon);
    treeTaggerInput.Copy(recon);
    treeTriggerInput.Copy(recon);
    treeTrackInput.Copy(recon);

    ++current_entry;
    return true;
}

double GoatReader::PercentDone() const
{
    return double(current_entry)/double(max_entries);
}

bool GoatReader::treeDetectorHitInput_t::LinkBranches(const WrapTFileInput& input, trees_t& trees)
{
    TTree* tree;
    if(!input.GetObject("detectorHits", tree))
        return false;
    NaI.LinkBranches(tree);
    PID.LinkBranches(tree);
    MWPC.LinkBranches(tree);
    BaF2.LinkBranches(tree);
    Veto.LinkBranches(tree);
    insert_trees(trees, NaI, PID, MWPC, BaF2, Veto);
    return true;
}

bool GoatReader::treeTaggerInput_t::LinkBranches(const WrapTFileInput& input, trees_t& trees)
{
    if(!input.GetObject("tagger",t.Tree))
        return false;
    t.LinkBranches();
    insert_trees(trees, t);
    return true;
}

bool GoatReader::treeTriggerInput_t::LinkBranches(const WrapTFileInput& input, trees_t& trees)
{
    if(!input.GetObject("trigger",t.Tree))
        return false;
    if(!input.GetObject("eventParameters",tEventParams.Tree))
        return false;
    t.LinkBranches();
    tEventParams.LinkBranches();

    if(t.MC_evt_id.IsPresent ^ t.MC_rnd_id.IsPresent)
        throw Exception("Branch MC_evt_id and MC_rnd_id inconsistenly present");

    LOG_IF(!t.helicity.IsPresent,  WARNING) << "Helicity bit information not found in input";
    LOG_IF(!t.MC_evt_id.IsPresent, WARNING) << "MC_evt_id/MC_rnd_id not found in input";

    insert_trees(trees, t, tEventParams);
    return true;
}

bool GoatReader::treeTrackInput_t::LinkBranches(const WrapTFileInput& input, trees_t& trees)
{
    if(!input.GetObject("tracks", t.Tree))
        return false;
    t.LinkBranches();
    insert_trees(trees, t);
    return true;
}



void GoatReader::treeDetectorHitInput_t::Copy(TEventData& recon)
{
    auto fill_readhits = [] (TEventData& recon,
            const treeEnergyTime_t& tree,
            Detector_t::Type_t type) {
        for(int i=0;i<int(tree.Hits().size());i++) {
            const unsigned channel = tree.Hits[i];
            const double energy = tree.Energy[i];
            const double time = tree.Time[i];
            if(isfinite(energy)) {
                recon.DetectorReadHits.emplace_back(
                            LogicalChannel_t{type, Channel_t::Type_t::Integral, channel},
                            TDetectorReadHit::Value_t{energy}
                            );
            }
            if(isfinite(time)) {
                recon.DetectorReadHits.emplace_back(
                            LogicalChannel_t{type, Channel_t::Type_t::Timing, channel},
                            TDetectorReadHit::Value_t{time}
                            );
            }
        }
    };

    fill_readhits(recon, NaI, Detector_t::Type_t::CB);
    fill_readhits(recon, PID, Detector_t::Type_t::PID);
    fill_readhits(recon, BaF2, Detector_t::Type_t::TAPS);
    fill_readhits(recon, Veto, Detector_t::Type_t::TAPSVeto);
    /// \todo think about MWPC stuff here, only hits are provided in GoAT tree?!
}



void GoatReader::treeTaggerInput_t::Copy(TEventData& recon)
{
    for( Int_t i=0; i<int(t.taggedChannel().size()); ++i) {
        recon.TaggerHits.emplace_back(
                    t.taggedChannel[i],
                    t.taggedEnergy[i],
                    t.taggedTime[i]
                    );
    }
}

void GoatReader::treeTriggerInput_t::Copy(TEventData& recon)
{
    TTrigger& ti = recon.Trigger;

    ti.DAQEventID = tEventParams.eventNumber;
    // we assume stuff from GoAT as measured, just to compare it with our
    // own simulation for this...
    ti.CBEnergySum = t.energySum; // not set as it's not really measured!
    ti.ClusterMultiplicity = t.multiplicity; // check if that's calculated as well?

    for( int err=0; err < int(t.errorCode().size()); ++err) {
        ti.DAQErrors.emplace_back(
                    t.errorModuleID[err],
                    t.errorModuleIndex[err],
                    t.errorCode[err]
                    );
    }
}

void GoatReader::treeTrackInput_t::Copy(TEventData& recon)
{
    /**
     * @brief map the cluster sizes from goat to unsigned ints
     * negative values mean no hit in the calorimeter
     * map those to 0
     */
    auto MapClusterSize =  [] (const int& size) {
        return size < 0 ? 0 : size;
    };

    /**
     * @brief map goat apparatus numbers to apparatus_t enum values
     * in case unknown values show up: -> exception and do not sliently ignore
     */
    auto IntToDetector_t =  [] (const int& a) {
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
    };


    for(Int_t i=0; i< int(t.clusterEnergy().size()); ++i) {

        const Detector_t::Any_t det = IntToDetector_t(t.detectors[i]);

        // Goat does not provide clusters,
        // so simulate some with fuzzy logic...
        TClusterList clusters;
        /// \todo how does this work with MWPC?

        if(det & Detector_t::Any_t::Calo) {

            clusters.emplace_back(
                        vec3::RThetaPhi(1.0, t.theta[i], t.phi[i]),
                        t.clusterEnergy[i],
                        t.time[i],
                        det & Detector_t::Type_t::CB ? Detector_t::Type_t::CB : Detector_t::Type_t::TAPS ,
                        t.centralCrystal[i]
                        );
            auto& calo_cluster = clusters.back();
            if(t.shortEnergy[i]>0)
                calo_cluster.ShortEnergy = t.shortEnergy[i];

            double vetoEnergy = 0.0;
            if(det & Detector_t::Any_t::Veto) {
                vetoEnergy =  t.vetoEnergy[i];
                clusters.emplace_back(
                            vec3(std_ext::NaN, std_ext::NaN, std_ext::NaN), // no veto position available
                            vetoEnergy,
                            std_ext::NaN, // no veto timing available
                            det & Detector_t::Type_t::PID ? Detector_t::Type_t::PID : Detector_t::Type_t::TAPSVeto,
                            t.centralVeto[i]
                            );

            }

            recon.Candidates.emplace_back(
                        det,
                        t.clusterEnergy[i],
                        t.theta[i],
                        t.phi[i],
                        t.time[i],
                        MapClusterSize(t.clusterSize[i]),
                        vetoEnergy,
                        //tracks.GetMWPC0Energy[i]+tracks.GetMWPC1Energy[i],
                        std_ext::NaN, // MWPC not handled correctly at the moment
                        TClusterList(clusters.begin(), clusters.end())
                        );
        }
        else if(det & Detector_t::Any_t::Veto) {
            // veto-only track is just a cluster in Ant
            // don't know if such tracks actually exist in GoAT/Acqu...
            const double vetoEnergy =  t.vetoEnergy[i];
            clusters.emplace_back(
                        vec3(std_ext::NaN, std_ext::NaN, std_ext::NaN), // no veto position available
                        vetoEnergy,
                        std_ext::NaN, // no veto timing available
                        det & Detector_t::Type_t::PID ? Detector_t::Type_t::PID : Detector_t::Type_t::TAPSVeto,
                        t.centralVeto[i]
                        );
        }

        // always add clusters...
        recon.Clusters.insert(recon.Clusters.end(), clusters.begin(), clusters.end());
    }
}
