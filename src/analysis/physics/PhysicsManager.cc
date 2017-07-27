#include "PhysicsManager.h"

#include "utils/ParticleID.h"
#include "input/DataReader.h"

#include "tree/TSlowControl.h"
#include "base/Logger.h"

#include "slowcontrol/SlowControlManager.h"

#include "base/ProgressCounter.h"

#include "TTree.h"

#include <iomanip>


using namespace std;
using namespace ant;
using namespace ant::analysis;

PhysicsManager::PhysicsManager(volatile bool* interrupt_) :
    physics(),
    interrupt(interrupt_),
    processedTIDrange(TID(), TID())
{}

PhysicsManager::~PhysicsManager() {}

void PhysicsManager::ShowResults()
{
    for(auto& p : physics) {
        p->ShowResult();
    }
}

void PhysicsManager::InitReaders(readers_t readers_)
{
    amenders = move(readers_);

    // figure out what set of readers we have
    // we expect maximum one source and several amenders
    auto it_amender = amenders.begin();
    while(it_amender != amenders.end()) {
        const auto flags = (*it_amender)->GetFlags();
        reader_flags |= flags;
        if(flags & input::reader_flag_t::IsSource) {
            if(source != nullptr)
                throw Exception("Found more than one source in given readers");
            source = move(*it_amender);
            it_amender = amenders.erase(it_amender);
        }
        else {
            ++it_amender;
        }
    }
}


void PhysicsManager::ReadFrom(
        list<unique_ptr<input::DataReader> > readers_,
        long long maxevents)
{
    InitReaders(move(readers_));

    if(physics.empty())
        throw Exception("No analysis instances activated. Cannot not analyse anything.");

    // prepare slowcontrol, init here since physics classes
    // register slowcontrol variables in constructor
    SlowControlManager slowControlManager(reader_flags);


    // prepare output of TEvents
    treeEvents.CreateBranches(new TTree("treeEvents","TEvent data"));

    long long nEventsRead = 0;
    long long nEventsProcessed = 0;
    long long nEventsAnalyzed = 0;
    long long nEventsSaved = 0;

    bool reached_maxevents = false;


    ProgressCounter progress(
                [this, &nEventsAnalyzed, maxevents]
                (std::chrono::duration<double> elapsed)
    {
        if (!source)
            return;
        const double percent = maxevents == numeric_limits<decltype(maxevents)>::max() ?
                                   source->PercentDone() :
                                   (double)nEventsAnalyzed/maxevents;

        static double last_PercentDone = 0;
        const double speed = (percent - last_PercentDone)/elapsed.count();
        LOG(INFO) << setw(2) << std::setprecision(4)
                  << percent*100 << " % done, ETA: " << ProgressCounter::TimeToStr((1-percent)/speed);
        last_PercentDone = percent;
    });
    while(true) {
        if(reached_maxevents || interrupt)
            break;

        // read events until slowcontrol_mgr is complete
        while(true) {
            if(interrupt) {
                VLOG(3) << "Reading interrupted";
                break;
            }

            input::event_t event;
            if(!TryReadEvent(event)) {
                VLOG(5) << "No more events to read, finish.";
                reached_maxevents = true;
                break;
            }
            nEventsRead++;

            // dump it into slowcontrol until full...
            if(slowControlManager.ProcessEvent(move(event)))
                break;
            // ..or max buffersize reached: 20000 corresponds to two Acqu Scaler blocks
            if(slowControlManager.BufferSize()>20000) {
                throw Exception(std_ext::formatter() <<
                                "Slowcontrol buffer reached maximum size " << slowControlManager.BufferSize()
                                << " without becoming complete. Stopping.");
            }
        }

        // read the slowcontrol_mgr's buffer and process the events
        while(auto buf_event = slowControlManager.PopEvent()) {

            auto& event = buf_event.Event;

            if(interrupt) {
                VLOG(3) << "Processing interrupted";
                break;
            }

            logger::DebugInfo::nProcessedEvents = nEventsProcessed;

            physics::manager_t manager;

            // if we've already reached the maxevents,
            // we just postprocess the remaining slowcontrol buffer (if any)
            if(!reached_maxevents) {
                if(nEventsAnalyzed == maxevents) {
                    VLOG(3) << "Reached max Events " << maxevents;
                    reached_maxevents = true;
                    // we cannot simply break here since might
                    // need to save stuff for slowcontrol purposes
                    if(slowControlManager.BufferSize()==0)
                        break;
                }

                if(!reached_maxevents && !buf_event.WantsSkip) {

                    ProcessEvent(event, manager);

                    // prefer Reconstructed ID, but at least one branch should be non-null
                    const auto& eventid = event.HasReconstructed() ? event.Reconstructed().ID : event.MCTrue().ID;
                    if(nEventsAnalyzed==0)
                        processedTIDrange.Start() = eventid;
                    processedTIDrange.Stop() = eventid;

                    nEventsAnalyzed++;

                    if(manager.saveEvent)
                        nEventsSaved++;
                }
            }

            // SaveEvent is the sink for events
            SaveEvent(move(event), manager);

            nEventsProcessed++;
        }
        ProgressCounter::Tick();
    }

    for(auto& pclass : physics) {
        pclass->Finish();
    }

    VLOG(5) << "Processed TID range: " << processedTIDrange;

    string processed_str;
    if(nEventsProcessed != nEventsAnalyzed)
        processed_str += std_ext::formatter() << " (" << nEventsProcessed << " processed)";
    if(nEventsRead != nEventsAnalyzed)
        processed_str += std_ext::formatter() << " (" << nEventsRead << " read)";


    LOG(INFO) << "Analyzed " << nEventsAnalyzed << " events"
              << processed_str << ", speed "
              << nEventsProcessed/progress.GetTotalSecs() << " event/s";

    const auto nEventsSavedTotal = treeEvents.Tree->GetEntries();
    if(nEventsSaved==0) {
        if(nEventsSavedTotal>0)
            VLOG(5) << "Deleting " << nEventsSavedTotal << " treeEvents from slowcontrol only";
        delete treeEvents.Tree;
    }
    else if(treeEvents.Tree->GetCurrentFile() != nullptr) {
        treeEvents.Tree->Write();
        const auto n_sc = nEventsSavedTotal - nEventsSaved;
        LOG(INFO) << "Wrote " << nEventsSaved  << " treeEvents"
                  << (n_sc>0 ? string(std_ext::formatter() << " (+slowcontrol: " << n_sc << ")") : "")
                  << ": "
                  << (double)treeEvents.Tree->GetTotBytes()/(1 << 20) << " MB (uncompressed), "
                  << (double)treeEvents.Tree->GetTotBytes()/nEventsSavedTotal << " bytes/event";
    }

    // cleanup readers (important for stopping progress output)
    source = nullptr;
    amenders.clear();
}


bool PhysicsManager::TryReadEvent(input::event_t& event)
{
    bool event_read = false;
    if(source) {
        if(!source->ReadNextEvent(event)) {
            return false;
        }
        event_read = true;
    }


    auto it_amender = amenders.begin();
    while(it_amender != amenders.end()) {

        if((*it_amender)->ReadNextEvent(event)) {
            ++it_amender;
            event_read = true;
        }
        else {
            it_amender = amenders.erase(it_amender);
        }
    }

    return event_read;
}

void PhysicsManager::ProcessEvent(input::event_t& event, physics::manager_t& manager)
{

    event.EnsureTempBranches();

    // run the physics classes
    for( auto& m : physics ) {
        m->ProcessEvent(event, manager);
    }

    event.ClearTempBranches();
}

void PhysicsManager::SaveEvent(input::event_t event, const physics::manager_t& manager)
{
    if(manager.saveEvent || event.SavedForSlowControls) {
        // only warn if manager says it should save
        if(!treeEvents.Tree->GetCurrentFile() && manager.saveEvent)
            LOG_N_TIMES(1, WARNING) << "Writing treeEvents to memory. Might be a lot of data!";


        // always keep read hits if saving for slowcontrol
        if(!manager.keepReadHits && !event.SavedForSlowControls)
            event.ClearDetectorReadHits();

        treeEvents.data = move(event);
        treeEvents.Tree->Fill();
    }
}
