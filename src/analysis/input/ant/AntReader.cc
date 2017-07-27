#include "AntReader.h"

#include "tree/TEvent.h"
#include "tree/TEventData.h"

#include "base/Logger.h"
#include "base/WrapTTree.h"
#include "input/treeEvents_t.h"

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
    virtual event_t NextEvent() = 0;
    virtual bool ProvidesSlowControl() const = 0;
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
    virtual event_t NextEvent() override {
        return event_t{unpacker->NextEvent()};
    }
    virtual bool ProvidesSlowControl() const override {
        return unpacker->ProvidesSlowControl();
    }
private:
    unique_ptr<Unpacker::Module> unpacker;
}; // UnpackerReader


struct TreeReader : AntReaderInternal {
    TreeReader(const std::shared_ptr<WrapTFileInput>& rootfiles)
    {
        if(!rootfiles->GetObject("treeEvents", tree.Tree))
            return;

        VLOG(5) << "Found Ant Events Tree";
        tree.LinkBranches();
    }

    virtual ~TreeReader() = default;

    virtual double PercentDone() const override {
        if(tree)
            return double(current_entry)/double(tree.Tree->GetEntries());
        return numeric_limits<double>::quiet_NaN();
    }

    virtual event_t NextEvent() override {
        if(!tree)
            return {};

        if(current_entry==tree.Tree->GetEntries())
            return {};

        tree.Tree->GetEntry(current_entry);
        current_entry++;
        return event_t{move(tree.data())};
    }

    virtual bool ProvidesSlowControl() const override {
        /// \todo the current implementation of reader flags and slow control providers looks non-optimal,
        /// improve this...
        LOG(WARNING) << "Reading from Ant trees assumes slow control present, this might fail if reading from MC-based trees.";
        return true;
    }

private:
    Long64_t current_entry = 0;

    treeEvents_t tree;
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

reader_flags_t AntReader::GetFlags() const {
    if(reader) {
        reader_flags_t flags(reader_flag_t::IsSource);
        if(reader->ProvidesSlowControl())
            flags |= reader_flag_t::ProvidesSlowControl;
        return flags;
    }
    else
        return {};
}

double AntReader::PercentDone() const
{
    if(reader)
        return reader->PercentDone();
    return numeric_limits<double>::quiet_NaN();
}

bool AntReader::ReadNextEvent(event_t& event)
{
    if(!reader)
        return false;

    // we expect Reconstructed branch to be filled always
    auto nextevent = reader->NextEvent();

    if(nextevent) {
        if(reconstruct) {
            TEventData& recon = nextevent.Reconstructed();
            /// \todo improve check if TEvent was run through reconstructed
            /// you may also introduce some flag to force application?
            if(recon.Clusters.empty())
                reconstruct->DoReconstruct(recon);
        }

        // pay attention that Geant unpacker might also set MCTrue branch partly
        event = move(nextevent);

        return true;
    }

    reader = nullptr;
    return false;
}

