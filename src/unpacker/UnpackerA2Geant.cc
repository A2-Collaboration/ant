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
            weights.emplace_back(
                        tagger->IsIgnored(ch) ? 0.0 : 1.0/tagger->GetPhotonEnergy(ch)
                        );
        }
        return {weights.begin(), weights.end()};
    }
};

}}}

using namespace ant::unpacker::geant;

UnpackerA2Geant::UnpackerA2Geant() {}

UnpackerA2Geant::~UnpackerA2Geant() {
    delete id;
}

bool UnpackerA2Geant::OpenFile(const string& filename)
{
    // open a root file, ignore error silently
    inputfile = std_ext::make_unique<WrapTFileInput>();

    try {
        inputfile->OpenFile(filename);
    } catch (const std::runtime_error&) {
        return false;
    }

    // setup the "expected" A2 geant tree
    if(!inputfile->GetObject("h12", geant))
        return false;

    geant->SetBranchAddress("nhits",&fnhits);
    geant->SetBranchAddress("npart",&fnpart);
    geant->SetBranchAddress("ntaps",&fntaps);
    geant->SetBranchAddress("nvtaps",&fnvtaps);
    geant->SetBranchAddress("vhits",&fvhits);
    geant->SetBranchAddress("plab",plab);
    geant->SetBranchAddress("tctaps",tctaps);
    geant->SetBranchAddress("vertex",fvertex);
    geant->SetBranchAddress("beam",fbeam);
    geant->SetBranchAddress("dircos",dircos);
    geant->SetBranchAddress("ecryst",ecryst);
    geant->SetBranchAddress("tcryst",tcryst);
    geant->SetBranchAddress("ectapfs",ectapfs);
    geant->SetBranchAddress("ectapsl",ectapsl);
    geant->SetBranchAddress("elab",elab);
    geant->SetBranchAddress("eleak",&feleak);
    geant->SetBranchAddress("enai",&fenai);
    geant->SetBranchAddress("etot",&fetot);
    geant->SetBranchAddress("eveto",eveto);
    geant->SetBranchAddress("tveto",tveto);
    geant->SetBranchAddress("evtaps",evtaps);
    geant->SetBranchAddress("icryst",icryst);
    geant->SetBranchAddress("ictaps",ictaps);
    geant->SetBranchAddress("ivtaps",ivtaps);
    geant->SetBranchAddress("idpart",idpart);
    geant->SetBranchAddress("iveto",iveto);
    geant->SetBranchAddress("nmwpc",&fnmwpc);
    geant->SetBranchAddress("imwpc",imwpc);
    geant->SetBranchAddress("mposx",mposx);
    geant->SetBranchAddress("mposy",mposy);
    geant->SetBranchAddress("mposz",mposz);
    geant->SetBranchAddress("emwpc",emwpc);

    TTree* tid_tree = nullptr;

    if(inputfile->GetObject("h12_tid", tid_tree)) {
        if(tid_tree->GetEntries() != geant->GetEntries()) {
            throw Exception("Geant Tree and TID Tree size mismatch");
        }

        geant->AddFriend(tid_tree);
        geant->SetBranchAddress("tid", &id);

        tid_from_file = true;

    } else {

        tid_from_file = false;

        /// \todo think of some better timestamp?
        id = new TID(static_cast<std::uint32_t>(std::time(nullptr)),
                     0, // start with 0 as lower ID
                     std::list<TID::Flags_t>{TID::Flags_t::MC, TID::Flags_t::AdHoc} // mark as MC
                     );
    }

    if(geant->GetEntries() >= numeric_limits<std::uint32_t>::max()) {
        throw Exception("Tree file contains too many entries for building correct unique ID");
    }

    geant->GetEntry(0);


    // try to get a config
    auto setup = ExpConfig::Setup::GetLastFound();
    if(!setup) {
        throw ExpConfig::ExceptionNoConfig("No setup found. This unpacker requires a manually set setup name");
    }
    auto config = dynamic_pointer_cast<UnpackerA2GeantConfig, ExpConfig::Setup>(setup);
    if(!config) {
        throw ExpConfig::ExceptionNoConfig("Found setup cannot configure this unpacker.");
    }


    // find some taggerdetectors
    // needed to create proper tagger hits from incoming photons
    for(const shared_ptr<Detector_t>& detector : config->GetDetectors()) {
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
                           config->GetPromptRandomConfig()
                           );
    }




    LOG(INFO) << "Successfully opened '" << filename
              << "' with " << geant->GetEntries() << " entries"
              << (tid_from_file ? ", with TID match check" : "");
    LOG_IF(!tid_from_file, WARNING) << "No TID match check enabled";

    return true;
}

struct r_t {
    std::default_random_engine r_gen;
    std::normal_distribution<double> r_prompt;
};


TEvent UnpackerA2Geant::NextEvent() noexcept
{
    if(current_entry>=geant->GetEntriesFast()-1)
        return {};

    geant->GetEntry(++current_entry);

    // start with an empty event with reconstructed ID set
    // MCTrue ID will be set by MCTrue reader, but this unpacker
    // knows the true vertex position...
    TEvent event(*id, TID());

    // however, vertex is some MCTrue information!
    event.MCTrue().Target.Vertex = vec3(fvertex[0], fvertex[1], fvertex[2]);

    const auto n_total = static_cast<unsigned>(fnhits+fnpart+fntaps+fnvtaps+fvhits);

    // approx. 3 detector read hits per detector, we just want to prevent re-allocation
    vector<TDetectorReadHit>& hits = event.Reconstructed().DetectorReadHits;
    hits.reserve(3*n_total);

    // all energies from A2geant are in GeV, but here we need MeV...
    const double GeVtoMeV = 1000.0;

    // fill CB Hits
    for(int i=0;i<fnhits;i++) {
        const auto ch = static_cast<unsigned>(icryst[i]); // no -1 here!

        if(ch >= cb_detector->GetNChannels())
            throw Exception("CB channel number out of bounds " + to_string(ch) + " / " + to_string(cb_detector->GetNChannels()));

        if(cb_detector->IsIgnored(ch))
            continue;
        const Detector_t::Type_t det = Detector_t::Type_t::CB;
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Integral, ch},
                    vector<double>{GeVtoMeV*ecryst[i]}
                    );
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Timing, ch},
                    vector<double>{tcryst[i]}
                    );
    }

    // fill PID Hits
    for(int i=0;i<fvhits;i++) {
        /// @todo Make PID channel mapping/rotation a Setup option?
        const unsigned ch = (23 - (iveto[i]-1) + 11) % 24;

        if(ch >= pid_detector->GetNChannels())
            throw Exception("PID channel number out of bounds " + to_string(ch) + " / " + to_string(pid_detector->GetNChannels()));

        if(pid_detector->IsIgnored(ch))
            continue;
        const Detector_t::Type_t det = Detector_t::Type_t::PID;
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Integral, ch},
                    vector<double>{GeVtoMeV*eveto[i]}
                    );
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Timing, ch},
                    vector<double>{tveto[i]}
                    );
    }

    // fill TAPS Hits
    for(int i=0;i<fntaps;i++) {
        const auto ch = static_cast<unsigned>(ictaps[i]-1);

        if(ch >= taps_detector->GetNChannels())
            throw Exception("TAPS channel number out of bounds " + to_string(ch) + " / " + to_string(taps_detector->GetNChannels()));

        if(taps_detector->IsIgnored(ch))
            continue;
        const Detector_t::Type_t det = Detector_t::Type_t::TAPS;
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Integral, ch},
                    vector<double>{GeVtoMeV*ectapsl[i]}
                    );
        /// \todo check if the short gate actually makes sense?
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::IntegralShort, ch},
                    vector<double>{GeVtoMeV*ectapfs[i]}
                    );
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Timing, ch},
                    vector<double>{tctaps[i]}
                    );
    }

    // fill TAPSVeto Hits
    for(int i=0;i<fnvtaps;i++) {
        const auto ch = static_cast<unsigned>(ivtaps[i]-1);

        if(ch >= tapsveto_detector->GetNChannels())
            throw Exception("TAPS channel number out of bounds " + to_string(ch) + " / " + to_string(tapsveto_detector->GetNChannels()));

        if(tapsveto_detector->IsIgnored(ch))
            continue;
        const Detector_t::Type_t det = Detector_t::Type_t::TAPSVeto;
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Integral, ch},
                    vector<double>{GeVtoMeV*evtaps[i]}
                    );
        /// \todo check if there's really no veto timing?
        hits.emplace_back(
                    LogicalChannel_t{det, Channel_t::Type_t::Timing, ch},
                    vector<double>{0}
                    );
    }

    // "reconstruct" a tagger electron from the photon
    const double photon_energy = GeVtoMeV*fbeam[4];

    if(taggerdetector) {
        // could the prompt photon have been detected?
        unsigned ch;
        if(taggerdetector->TryGetChannelFromPhoton(photon_energy, ch) &&
           !taggerdetector->IsIgnored(ch)
           )
        {
            // then insert (possibly time-smeared) prompt hit
            hits.emplace_back(
                        LogicalChannel_t{taggerdetector->Type, Channel_t::Type_t::Timing, ch},
                        std::vector<double>{promptrandom->SmearPrompt(0)}
                        );


        }

        // always fill some extra random hits
        for(auto& hit : promptrandom->GetRandomHits()) {
            hits.emplace_back(
                        LogicalChannel_t{taggerdetector->Type, Channel_t::Type_t::Timing, hit.Channel},
                        std::vector<double>{hit.Timing}
                        );
        }
    }


    if(!tid_from_file)
        ++(*id);

    return event;
}

double UnpackerA2Geant::PercentDone() const
{
    return double(current_entry) / double(geant->GetEntries());
}


