#include "AntReader.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "base/Logger.h"

#include "TTree.h"

#include <memory>
#include <stdexcept>

using namespace std;
using namespace ant;
using namespace ant::analysis::input;

namespace ant {
namespace analysis {
namespace input {
namespace detail {

struct AntReaderInternal {
    virtual double PercentDone() const = 0;
    virtual std::unique_ptr<TEvent> NextEvent() = 0;
    virtual ~AntReaderInternal() = default;
};


struct UnpackerReader : AntReaderInternal {
    UnpackerReader(unique_ptr<Unpacker::Module> unpacker_) :
        unpacker(move(unpacker_))
    {
        LOG(INFO) << "Reading events from unpacker";
    }
    virtual double PercentDone() const override {
        return unpacker->PercentDone();
    }
    virtual std::unique_ptr<TEvent> NextEvent() override {
        return unpacker->NextEvent();
    }
private:
    unique_ptr<Unpacker::Module> unpacker;
}; // UnpackerReader


struct TreeReader : AntReaderInternal {
    TreeReader(const std::shared_ptr<WrapTFileInput>& rootfiles)
    {
        if(!rootfiles->GetObject("treeEvents", tree))
            return;

        VLOG(5) << "Found Ant Events Tree";

        const auto res = tree->SetBranchAddress("data", addressof(eventPtr));
        if(res != TTree::kMatch) {
            tree = nullptr;
            LOG(ERROR) << "Could not access branch 'data' in Ant events tree";
            return;
        }
    }

    virtual double PercentDone() const override {
        if(tree)
            return double(current_entry)/double(tree->GetEntries());
        return numeric_limits<double>::quiet_NaN();
    }

    virtual std::unique_ptr<TEvent> NextEvent() override {
        if(!tree)
            return nullptr;

        if(current_entry==tree->GetEntries())
            return nullptr;

        eventPtr = nullptr;
        tree->GetEntry(current_entry);
        current_entry++;
        return unique_ptr<TEvent>(eventPtr);
    }

private:
    Long64_t current_entry = 0;
    TTree* tree = nullptr;
    TEvent* eventPtr = nullptr;

}; // TreeReader

}}}} // namespace ant::analysis::input::detail


AntReader::AntReader(const std::shared_ptr<WrapTFileInput>& rootfiles,
        unique_ptr<Unpacker::Module> unpacker,
        std::unique_ptr<Reconstruct_traits> reconstruct_
        ) :
    reconstruct(move(reconstruct_))
{
    // prefer unpacker
    if(unpacker) {
        reader = std_ext::make_unique<detail::UnpackerReader>(move(unpacker));
        if(!reconstruct)
            LOG(WARNING) << "Reconstruct disabled although reading from unpacker. Producing DetectorReadHits only.";
    }
    else {
        // try root files
        auto treereader = std_ext::make_unique<detail::TreeReader>(rootfiles);
        if(isfinite(treereader->PercentDone()))
            reader = move(treereader);
    }

}

AntReader::~AntReader() {}

bool AntReader::IsSource() {
    return reader != nullptr;
}

double AntReader::PercentDone() const
{
    if(reader)
        return reader->PercentDone();
    return numeric_limits<double>::quiet_NaN();
}

bool AntReader::ReadNextEvent(TEvent& event)
{
    if(!reader)
        return false;

    // we expect Reconstructed branch to be filled always
    auto eventptr = reader->NextEvent();

    if(eventptr) {
        if(reconstruct) {
            TEventData& recon = eventptr->Reconstructed();
            /// \todo improve check if TEvent was run through reconstructed
            /// you may also introduce some flag to force application?
            if(recon.Clusters.empty())
                reconstruct->DoReconstruct(recon);
        }

        // pay attention that Geant unpacker might also set MCTrue branch partly
        event = move(*eventptr);

        return true;
    }

    reader = nullptr;
    return false;
}

