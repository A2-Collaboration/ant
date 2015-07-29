

#include "analysis/input/DataReader.h"
#include "analysis/input/ant/AntReader.h"
#include "analysis/input/goat/GoatReader.h"
#include "analysis/OutputManager.h"

#include "analysis/physics/Physics.h"
#include "analysis/physics/omega/omega.h"
#include "analysis/physics/common/DataOverview.h"
#include "analysis/physics/common/CandidatesAnalysis.h"

#include "expconfig/ExpConfig.h"

#include "unpacker/Unpacker.h"

#include "tree/THeaderInfo.h"
#include "tree/TUnpackerMessage.h"
#include "tree/TSlowControl.h"
#include "tree/TDetectorRead.h"
#include "tree/TEvent.h"

#include "base/std_ext.h"
#include "base/Logger.h"
#include "base/CmdLine.h"
#include "base/ReadTFiles.h"

#include "TRint.h"
#include "TClass.h"
#include "TTree.h"
#include "TFile.h"

#include <sstream>
#include <string>

using namespace std;
using namespace ant::output;
using namespace ant;
using namespace ant::analysis;



class UnpackerTreeReader {

    shared_ptr<ReadTFiles> file;

    TEvent* Event;
    TDetectorRead* DetectorRead;
    THeaderInfo* HeaderInfo;
    TUnpackerMessage* UnpackerMessage;
    TSlowControl* SlowControl;

    struct treerecord_t {
        TDataRecord** Record;
        Long64_t CurrEntry;
        TTree* Tree;
        const TDataRecord& GetRecord() const {
            return *(*Record);
        }
    };

    using treerecords_t = vector<treerecord_t>;
    treerecords_t treerecords;
    treerecords_t::iterator it_treerecord;
    bool isopen;

    TID currID;

    TID findMinID() const {
        TID tid_min = it_treerecord->GetRecord().ID; // start with something existing
        for(const treerecord_t& treerecord : treerecords) {
            const TDataRecord& record = treerecord.GetRecord();
            if(record.ID < tid_min)
                tid_min = record.ID;
        }
        return tid_min;
    }

public:

    UnpackerTreeReader(const shared_ptr<ReadTFiles>& rootfiles) :
          file(rootfiles), // remember the shared_ptr to make sure it lives as long as this class
          treerecords(),
          it_treerecord(treerecords.end()), // invalid iterator position
          isopen(false),
          currID()
    {}

    template<typename T>
    bool SetupBranch(const string& name, T*& ptr) {
        TTree* tree = nullptr;
        const string treename = string("tree") + name;
        if(!file->GetObject(treename, tree))
            return false;
        if(tree->GetEntries()==0)
            return true;
        if(!tree->GetListOfBranches()->FindObject(name.c_str()))
            return false;
        tree->SetBranchAddress(name.c_str(), addressof(ptr));
        treerecords.emplace_back(treerecord_t{
                    reinterpret_cast<TDataRecord**>(addressof(ptr)),
                    0,
                    tree}
                    );
        ptr = new T();
        tree->GetEntry(0);
        return true;
    }

    bool OpenInput() {
        if(!SetupBranch("HeaderInfo", HeaderInfo))
            return false;
        if(!SetupBranch("UnpackerMessage", UnpackerMessage))
            return false;
        if(!SetupBranch("SlowControl", SlowControl))
            return false;
        if(!SetupBranch("Event", Event))
            return false;
        if(!SetupBranch("DetectorRead", DetectorRead))
            return false;

        if(treerecords.empty())
            return false;

        // set it_treerecord to the item with the minimum ID
        it_treerecord = treerecords.begin();
        currID = findMinID();
        while(it_treerecord != treerecords.end()) {
            if(currID == it_treerecord->GetRecord().ID)
                break;
            ++it_treerecord;
        }
        if(it_treerecord == treerecords.end())
            return false;

        VLOG(9) << "Start reading at ID=" << currID << " with type="
                << it_treerecord->GetRecord().IsA()->GetName();

        isopen = true;
        return true;
    }

    bool IsOpen() const {
        return isopen;
    }

    bool GetUniqueHeaderInfo(THeaderInfo& headerInfo) {
        for(const treerecord_t& treerecord : treerecords) {
            const TDataRecord& rec = treerecord.GetRecord();
            if(rec.IsA() == THeaderInfo::Class()) {
                TTree* tree = treerecord.Tree;
                if(tree->GetEntries() != 1)
                    return false;
                tree->GetEntry(0);
                headerInfo = *HeaderInfo;
                return true;
            }
        }
        return false;
    }

    virtual std::shared_ptr<TDataRecord> NextItem() noexcept {
        // the TID determines, what item is next,
        // if there are multiple, the current iterator decides


        return nullptr;

    }



//    void Fill(TDataRecord* rec) {
//        const TClass* isA = rec->IsA();
//        if(isA == TEvent::Class()) {
//            Event = static_cast<TEvent*>(rec);
//            treeEvent->Fill();
//        }
//        else if(isA == TDetectorRead::Class()) {
//            DetectorRead = static_cast<TDetectorRead*>(rec);
//            treeDetectorRead->Fill();
//        }
//        else if(isA == THeaderInfo::Class()) {
//            HeaderInfo = static_cast<THeaderInfo*>(rec);
//            treeHeaderInfo->Fill();
//        }
//    }



};

int main(int argc, char** argv) {
    SetupLogger();

    TCLAP::CmdLine cmd("ant", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"int");
    auto cmd_input  = cmd.add<TCLAP::MultiArg<string>>("i","input","Input files",true,"string");
    auto cmd_setup  = cmd.add<TCLAP::ValueArg<string>>("s","setup","Choose setup",false,"","string");



    cmd.parse(argc, argv);
    if(cmd_verbose->isSet()) {
        el::Loggers::setVerboseLevel(cmd_verbose->getValue());
    }

    // build the general file manager first
    auto filemanager = make_shared<ReadTFiles>();
    for(const auto& inputfile : cmd_input->getValue()) {
        VLOG(5) << "ROOT File Manager: Looking at file " << inputfile;
        if(filemanager->OpenFile(inputfile))
            LOG(INFO) << "Opened file '" << inputfile << "' as ROOT file";
    }

    // then init the unpacker root input file manager
    UnpackerTreeReader unpackerFile(filemanager);

    // search for header info?
    if(unpackerFile.OpenInput()) {
        LOG(INFO) << "Found complete set of input ROOT trees for unpacker";
        THeaderInfo headerInfo;
        if(unpackerFile.GetUniqueHeaderInfo(headerInfo)) {
            VLOG(5) << "Found unique header info " << headerInfo;
            if(!headerInfo.SetupName.empty()) {
                ExpConfig::ManualSetupName = headerInfo.SetupName;
                LOG(INFO) << "Using header info to manually set the setup name to " << ExpConfig::ManualSetupName;
            }
        }
    }

    // override the setup name from cmd line
    if(cmd_setup->isSet()) {
        ExpConfig::ManualSetupName = cmd_setup->getValue();
        if(ExpConfig::ManualSetupName.empty())
            LOG(INFO) << "Commandline override to auto-search for setup config (might fail)";
        else
            LOG(INFO) << "Commandline override setup name to '" << ExpConfig::ManualSetupName << "'";
    }


    // now we can try to open the files with the unpackers
    std::unique_ptr<Unpacker::Module> unpacker = nullptr;
    for(const auto& inputfile : cmd_input->getValue()) {
        VLOG(5) << "Unpacker: Looking at file " << inputfile;
        try {
            auto unpacker_ = Unpacker::Get(inputfile);
            if(unpacker != nullptr && unpacker_ != nullptr) {
                LOG(ERROR) << "Can only handle one unpacker, but input files suggest to use more than one.";
                return 1;
            }
            unpacker = move(unpacker_);
        }
        catch(Unpacker::Exception e) {
            VLOG(5) << "Unpacker::Get said: " << e.what();
            unpacker = nullptr;
        }
    }


    if(unpacker != nullptr && unpackerFile.IsOpen()) {
        LOG(WARNING) << "Found file suitable for unpacker and ROOT file for unpacker stage, preferring raw data file!";
    }


    //unpackerFile.DetectorRead->Class()->

    //cout << unpackerFile.HeaderInfo->Class_Name() << endl;

    return 0;

}

