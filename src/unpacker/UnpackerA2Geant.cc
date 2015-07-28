#include "UnpackerA2Geant.h"

#include "expconfig/ExpConfig.h"

#include "base/Logger.h"
#include "base/std_ext.h"

#include "tree/THeaderInfo.h"
#include "tree/TDetectorRead.h"

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

#include <memory>

using namespace std;
using namespace ant;

UnpackerA2Geant::UnpackerA2Geant() {}

UnpackerA2Geant::~UnpackerA2Geant() {
    tfile->Close();
}

bool UnpackerA2Geant::OpenFile(const string& filename)
{
    // open a root file, ignore error silently
    const auto prev_gErrorIgnoreLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError+1;
    tfile = std_ext::make_unique<TFile>(filename.c_str(),"READ");
    gErrorIgnoreLevel = prev_gErrorIgnoreLevel;
    if(tfile->IsZombie())
        return false;

    // setup the "expected" A2 geant tree
    geant = dynamic_cast<TTree*>(tfile->Get("h12"));
    if(geant == nullptr)
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
    geant->SetBranchAddress("mc_evt_id",&mc_evt_id);
    geant->SetBranchAddress("mc_rnd_id",&mc_evt_id);

    if(geant->GetEntries() >= numeric_limits<decltype(ID_lower)>::max()) {
        throw Exception("Tree file contains too many entries for building correct unique ID");
    }



    /// \todo think of some better upper Id?
    ID_lower = 0; // also counts number of entries in TTree
    ID_upper = geant->Hash();

    // this unpacker has no chance to make a proper THeaderInfo
    // so we ask the ExpConfig if it has an idea...

    if(ExpConfig::ManualSetupName.empty()) {
        throw Exception("This unpacker requires a manually set setup name");
    }
    headerInfo = std_ext::make_unique<THeaderInfo>(TID(ID_upper, ID_lower, true),
                                                   ExpConfig::ManualSetupName);
    config = ExpConfig::Unpacker<UnpackerA2GeantConfig>::Get(*headerInfo);

    LOG(INFO) << "Successfully opened tree in '" << filename
              << "' with " << geant->GetEntries() << " entries";

    return true;
}

shared_ptr<TDataRecord> UnpackerA2Geant::NextItem() noexcept
{
    if(ID_lower>=geant->GetEntriesFast())
        return nullptr;

    // return the headerinfo as the very first item
    if(headerInfo != nullptr) {
        return make_shared<THeaderInfo>(*headerInfo.release());
    }


    geant->GetEntry(ID_lower);

    // start with an empty detector read, don't forget MC flag true
    shared_ptr<TDetectorRead> detread = make_shared<TDetectorRead>(TID(ID_upper, ID_lower, true));

    const size_t n_total = fnhits+fnpart+fntaps+fnvtaps+fvhits;


    vector<TDetectorReadHit>& hits = detread->Hits;
    hits.reserve(3*n_total); // approx. 3 detector read hits per detector

    const double GeVtoMeV = 1000.0;

    // fill CB Hits
    for(int i=0;i<fnhits;i++) {
        const unsigned ch = icryst[i]; // no -1 here!
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
        /// \todo take care of reversed PID orientation
        const unsigned ch = iveto[i]-1;
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
        const unsigned ch = ictaps[i]-1;
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
        /// \todo check if -1 here is really correct
        const unsigned ch = ivtaps[i]-1;
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

    /// \todo simulate tagger hits...

    ID_lower++;
    return detread;
}


