#include "UnpackerA2Geant.h"

#include "expconfig/ExpConfig.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "base/WrapTFile.h"
#include "base/Logger.h"

#include "TTree.h"

#include <memory>
#include <random>

using namespace std;
using namespace ant;

namespace ant {
namespace unpacker {
namespace geant {

struct promptrandom_t {

    using tagger_t = shared_ptr<TaggerDetector_t>;

    promptrandom_t(
            tagger_t tagger,
            UnpackerA2GeantConfig::promptrandom_config_t config) :
        r_prompt(config.PromptOffset, config.PromptSigma),
        n_randoms(static_cast<unsigned>( // round-off error negligble
                      config.TimeWindow.Length()*config.RandomPromptRatio)
                  ),
        r_random_timing(config.TimeWindow.Start(), config.TimeWindow.Stop()),
        r_random_channel(make_r_tagg_ch(tagger))
    {
    }

    double SmearPrompt(double in) {
        return in + r_prompt(r_gen);
    }

    struct hit_t {
        unsigned Channel;
        double Timing;
    };

    std::vector<hit_t> GetRandomHits() {
        vector<hit_t> hits(n_randoms);
        for(unsigned i=0;i<n_randoms;i++) {
            auto& hit = hits[i];
            hit.Timing = r_random_timing(r_gen);
            hit.Channel = r_random_channel(r_gen);
        }
        return hits;
    }


protected:


    default_random_engine r_gen;

    normal_distribution<double> r_prompt;
    const unsigned n_randoms;
    uniform_real_distribution<double> r_random_timing;
    discrete_distribution<unsigned> r_random_channel;

    static decltype(r_random_channel) make_r_tagg_ch(tagger_t tagger) {
        /// \todo check if that channel weighting actually works here
        vector<double> weights;
        for(unsigned ch=0;ch<tagger->GetNChannels();ch++) {
            weights.emplace_back(1.0/tagger->GetPhotonEnergy(ch));
        }
        return {weights.begin(), weights.end()};
    }
};

}}}

using namespace ant::unpacker::geant;

UnpackerA2Geant::UnpackerA2Geant() {}

UnpackerA2Geant::~UnpackerA2Geant() {}

bool UnpackerA2Geant::OpenFile(const string& filename)
{
    // open a root file, ignore non-ROOT files silently
    inputfile = std_ext::make_unique<WrapTFileInput>();

    try {
        inputfile->OpenFile(filename);
    } catch (WrapTFile::ENotARootFile&) {
        return false;
    }

    // setup the "expected" A2 geant tree
    if(!inputfile->GetObject("h12", geantTree.Tree))
        return false;

    geantTree.LinkBranches();

    if(inputfile->GetObject("h12_tid", tidTree.Tree)) {
        if(tidTree.Tree->GetEntries() != geantTree.Tree->GetEntries()) {
            throw Exception("Geant Tree and TID Tree size mismatch");
        }
        tidTree.LinkBranches();
    } else {
        // think of some better timestamp?
        tidTree.tid = TID(static_cast<std::uint32_t>(std::time(nullptr)),
                          0, // start with 0 as lower ID
                          std::list<TID::Flags_t>{TID::Flags_t::MC, TID::Flags_t::AdHoc} // mark as MC
                          );
    }

    if(geantTree.Tree->GetEntries() >= numeric_limits<std::uint32_t>::max()) {
        throw Exception("Tree file contains too many entries for building proper unique ID");
    }

    // heuristically detect some older format, flag is used in unpacking
    {
        auto& t = geantTree;
        if(!t.tcryst.IsPresent && !t.tveto.IsPresent &&
           !t.ivtaps.IsPresent && !t.imwpc.IsPresent &&
           !t.mposx.IsPresent &&  !t.mposy.IsPresent &&
           !t.mposz.IsPresent && !t.emwpc.IsPresent) {
            LOG(INFO) << "Detected old format tree input, as some branches are missing.";
            oldTreeFormat = true;
        }
    }

    // try to get a config
    auto& setup = ExpConfig::Setup::GetByType<UnpackerA2GeantConfig>();

    // find some taggerdetectors
    // needed to create proper tagger hits from incoming photons
    for(const shared_ptr<Detector_t>& detector : setup.GetDetectors()) {
        /// \todo check for multiply defined detectors...
        if(auto tagger = dynamic_pointer_cast<TaggerDetector_t, Detector_t>(detector))
            taggerdetector = tagger;
        if(detector->Type == Detector_t::Type_t::CB)
            cb_detector = detector;
        if(detector->Type == Detector_t::Type_t::PID)
            pid_detector = detector;
        if(detector->Type == Detector_t::Type_t::TAPS)
            taps_detector = detector;
        if(detector->Type == Detector_t::Type_t::TAPSVeto)
            tapsveto_detector = detector;
    }

    if(!taggerdetector)
        LOG(WARNING) << "No tagger detector found in config, there will be no taggerhits generated";
    else {
        // initialize randomness
        promptrandom = std_ext::make_unique<unpacker::geant::promptrandom_t>(
                           taggerdetector,
                           setup.GetPromptRandomConfig()
                           );
    }




    LOG(INFO) << "Successfully opened '" << filename
              << "' with " << geantTree.Tree->GetEntries() << " entries"
              << (tidTree ? ", with TID match check" : "");
    LOG_IF(!tidTree, WARNING) << "No TID match check enabled";

    return true;
}

struct r_t {
    std::default_random_engine r_gen;
    std::normal_distribution<double> r_prompt;
};


TEvent UnpackerA2Geant::NextEvent()
{
    // shortcut, as geantTree is used very often here
    auto& t = geantTree;

    if(current_entry>=t.Tree->GetEntriesFast()-1)
        return {};

    t.Tree->GetEntry(++current_entry);

    // read TIDs in sync
    if(tidTree)
        tidTree.Tree->GetEntry(current_entry);

    // start with an empty event with reconstructed ID set
    // MCTrue ID will be set by MCTrue reader, but this unpacker
    // knows the true vertex position...
    TEvent event(tidTree.tid, TID());

    // however, vertex is some MCTrue information!
    event.MCTrue().Target.Vertex = vec3(t.vertex[0], t.vertex[1], t.vertex[2]);

    auto& hits = event.Reconstructed().DetectorReadHits;

    // all energies from A2geant are in GeV, but here we need MeV...
    const double GeVtoMeV = 1000.0;

    // fill CB Hits
    for(int i=0;i<int(t.icryst().size());i++) {
        const auto nCh = cb_detector->GetNChannels();
        if(oldTreeFormat) {
            if(t.icryst[i]<0 || t.icryst[i]>=static_cast<int>(nCh)) {
                LOG_N_TIMES(10, WARNING) << "Ignoring CB index out of bounds: " << t.icryst[i]
                                            << " i=" << i << " (max 10 times reported)";
                continue;
            }
        }

        const auto ch = static_cast<unsigned>(t.icryst[i]); // no -1 here!

        if(ch >= nCh)
            throw Exception("CB channel number out of bounds " + to_string(ch) + " / " + to_string(cb_detector->GetNChannels()));

        const Detector_t::Type_t det = Detector_t::Type_t::CB;
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Integral, ch},
                    TDetectorReadHit::Value_t{GeVtoMeV*t.ecryst[i]}
                    );
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Timing, ch},
                    TDetectorReadHit::Value_t{t.tcryst.IsPresent ? t.tcryst[i] : 0.0}
                    );
    }

    // fill PID Hits
    for(int i=0;i<int(t.iveto().size());i++) {
        /// @todo Make PID channel mapping/rotation a Setup option?
        const unsigned ch = (23 - (t.iveto[i]-1) + 11) % 24;

        if(ch >= pid_detector->GetNChannels())
            throw Exception("PID channel number out of bounds " + to_string(ch) + " / " + to_string(pid_detector->GetNChannels()));

        const Detector_t::Type_t det = Detector_t::Type_t::PID;
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Integral, ch},
                    TDetectorReadHit::Value_t{GeVtoMeV*t.eveto[i]}
                    );
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Timing, ch},
                    TDetectorReadHit::Value_t{t.tveto.IsPresent ? t.tveto[i] : 0.0}
                    );
    }

    // fill TAPS Hits
    for(int i=0;i<int(t.ictaps().size());i++) {
        // the older format appears to have some more "sane" index handling...
        const auto ch = static_cast<unsigned>(t.ictaps[i] - (oldTreeFormat ? 0 : 1));

        if(ch >= taps_detector->GetNChannels())
            throw Exception("TAPS channel number out of bounds " + to_string(ch) + " / " + to_string(taps_detector->GetNChannels()));

        const Detector_t::Type_t det = Detector_t::Type_t::TAPS;
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Integral, ch},
                    TDetectorReadHit::Value_t{GeVtoMeV*t.ectapsl[i]}
                    );
        /// \todo check if the short gate actually makes sense?
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::IntegralShort, ch},
                    TDetectorReadHit::Value_t{GeVtoMeV*t.ectapfs[i]}
                    );
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Timing, ch},
                    TDetectorReadHit::Value_t{t.tctaps[i]}
                    );
    }

    // fill TAPSVeto Hits
    if(t.ivtaps.IsPresent) {
        for(int i=0;i<int(t.ivtaps().size());i++) {
            const auto ch = static_cast<unsigned>(t.ivtaps[i]-1);

            if(ch >= tapsveto_detector->GetNChannels())
                throw Exception("TAPS channel number out of bounds " + to_string(ch) + " / " + to_string(tapsveto_detector->GetNChannels()));

            const Detector_t::Type_t det = Detector_t::Type_t::TAPSVeto;
            hits.emplace_back(
                        LogicalChannel_t{det, Channel_t::Type_t::Integral, ch},
                        TDetectorReadHit::Value_t{GeVtoMeV*t.evtaps[i]}
                        );
            /// \todo check if there's really no veto timing?
            hits.emplace_back(
                        LogicalChannel_t{det, Channel_t::Type_t::Timing, ch},
                        TDetectorReadHit::Value_t{0}
                        );
        }
    }

    // "reconstruct" a tagger electron from the photon
    const double photon_energy = GeVtoMeV*t.beam[4];

    if(taggerdetector) {
        // could the prompt photon have been detected?
        unsigned ch;
        if(taggerdetector->TryGetChannelFromPhoton(photon_energy, ch))
        {
            // then insert (possibly time-smeared) prompt hit
            hits.emplace_back(
                        LogicalChannel_t{taggerdetector->Type, Channel_t::Type_t::Timing, ch},
                        TDetectorReadHit::Value_t{promptrandom->SmearPrompt(0)}
                        );


        }

        // always fill some extra random hits
        for(auto& hit : promptrandom->GetRandomHits()) {
            hits.emplace_back(
                        LogicalChannel_t{taggerdetector->Type, Channel_t::Type_t::Timing, hit.Channel},
                        TDetectorReadHit::Value_t{hit.Timing}
                        );
        }
    }

    if(!tidTree)
        ++tidTree.tid();

    return event;
}

double UnpackerA2Geant::PercentDone() const
{
    return double(current_entry) / double(geantTree.Tree->GetEntries());
}


