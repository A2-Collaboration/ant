#include "PhysicsManager.h"

#include "utils/ParticleID.h"
#include "input/DataReader.h"

#include "tree/TSlowControl.h"
#include "tree/TAntHeader.h"
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
    interrupt(interrupt_)
{}

PhysicsManager::~PhysicsManager() {}

void PhysicsManager::SetParticleID(std::unique_ptr<utils::ParticleID> pid) {
    particleID = move(pid);
}

void PhysicsManager::SetAntHeader(TAntHeader& header)
{
    header.FirstID = firstID;
    header.LastID = lastID;
}

void PhysicsManager::EnableProgressUpdates(bool updates)
{
    progressUpdates = updates;
}

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
        if((*it_amender)->IsSource()) {
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
    slowcontrol_mgr = std_ext::make_unique<SlowControlManager>();


    // prepare output of TEvents
    treeEvents = new TTree("treeEvents","TEvent data");
    treeEventPtr = nullptr;
    treeEvents->Branch("data", addressof(treeEventPtr));


    chrono::time_point<std::chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    long long nEventsProcessed = 0;
    long long nEventsAnalyzed = 0;
    long long nEventsSaved = 0;

    bool reached_maxevents = false;

    ProgressCounter progress;
    while(true) {
        if(reached_maxevents || interrupt)
            break;

        // read events until slowcontrol_mgr is complete
        while(true) {
            if(interrupt) {
                VLOG(3) << "Reading interrupted";
                break;
            }

            auto event = std_ext::make_unique<TEvent>();
            if(!TryReadEvent(event)) {
                VLOG(5) << "No more events to read, finish.";
                reached_maxevents = true;
                break;
            }
            // dump it into slowcontrol until full...
            if(slowcontrol_mgr->ProcessEvent(move(event)))
                break;
            // ..or max buffersize reached: 20000 corresponds to two Acqu Scaler blocks
            if(slowcontrol_mgr->BufferSize()>20000) {
                throw Exception(std_ext::formatter() <<
                                "Slowcontrol buffer reached maximum size " << slowcontrol_mgr->BufferSize()
                                << " without becoming complete. Stopping.");
            }
        }

        // read the slowcontrol_mgr's buffer and process the events
        while(auto buffered_event = slowcontrol_mgr->PopEvent()) {

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
                    if(slowcontrol_mgr->BufferSize()==0)
                        break;
                }

                TEvent& event = *buffered_event.Event;
                if(!reached_maxevents && !event.SavedForSlowControls) {
                    ProcessEvent(event, manager);

                    // prefer Reconstructed ID, but at least one branch should be non-null
                    const auto& eventid = event.Reconstructed ? event.Reconstructed->ID : event.MCTrue->ID;
                    if(nEventsAnalyzed==0)
                        firstID = eventid;
                    lastID = eventid;

                    nEventsAnalyzed++;

                    if(manager.saveEvent)
                        nEventsSaved++;
                }
            }

            SaveEvent(move(buffered_event), manager);

            nEventsProcessed++;
        }

        if(progressUpdates && source && progress.Update(source->PercentDone())) {
            LOG(INFO) << progress;
        }

    }

    for(auto& pclass : physics) {
        pclass->Finish();
    }

    VLOG(5) << "First EventId processed: " << firstID;
    VLOG(5) << "Last  EventId processed: " << lastID;



    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    LOG(INFO) << "Analyzed " << nEventsAnalyzed << " events"
              << (
                     nEventsProcessed != nEventsAnalyzed ?
                                             string(std_ext::formatter() << " (" << nEventsProcessed << " processed)")
                                           : ""
                 )
              <<   ", speed "
              << nEventsProcessed/elapsed_seconds.count() << " event/s";

    const auto nEventsSavedTotal = treeEvents->GetEntries();
    if(nEventsSaved==0) {
        if(nEventsSavedTotal>0)
            VLOG(5) << "Deleting " << nEventsSavedTotal << " treeEvents from slowcontrol only";
        delete treeEvents;
    }
    else if(treeEvents->GetCurrentFile() != nullptr) {
        treeEvents->Write();
        const auto n_sc = nEventsSavedTotal - nEventsSaved;
        LOG(INFO) << "Wrote " << nEventsSaved  << " treeEvents"
                  << (n_sc>0 ? string(std_ext::formatter() << " (+slowcontrol: " << n_sc << ")") : "")
                  << ": "
                  << (double)treeEvents->GetTotBytes()/(1 << 20) << " MB (uncompressed), "
                  << (double)treeEvents->GetTotBytes()/nEventsSavedTotal << " bytes/event";
    }
}


bool PhysicsManager::TryReadEvent(TEventPtr& event)
{
    bool event_read = false;
    if(source) {
        if(!source->ReadNextEvent(*event)) {
            return false;
        }
        event_read = true;
    }


    auto it_amender = amenders.begin();
    while(it_amender != amenders.end()) {

        if((*it_amender)->ReadNextEvent(*event)) {
            ++it_amender;
            event_read = true;
        }
        else {
            it_amender = amenders.erase(it_amender);
        }
    }

    return event_read;
}

void PhysicsManager::ProcessEvent(TEvent& event, physics::manager_t& manager)
{
    if(particleID && event.Reconstructed) {
        // run particle ID for Reconstructed candidates
        // but only if there are no identified particles present yet
        /// \todo implement flag to force particle ID again?
        TEventData& recon = *event.Reconstructed;
        if(recon.Particles.GetAll().empty()) {
            for(const auto& cand : recon.Candidates) {
                auto particle = particleID->Process(cand);
                if(particle)
                    recon.Particles.Add(particle);
            }
        }
    }

    // ensure that physics classes always
    // have at least empty TEventData branches MCTrue and Reconstructed
    // use RAII for that
    struct clean_branch_t {
        TEventDataPtr& branch;
        bool clean = false;
        clean_branch_t(TEventDataPtr& branch_) : branch(branch_) {
            if(!branch) {
                // create temporary empty branch
                branch = std_ext::make_unique<TEventData>();
                clean = true;
            }
        }
        ~clean_branch_t() {
            // delete the branch if it was just empty in dtor
            if(clean)
                branch = nullptr;
        }
    };

    // make sure the branches are non-null for physics classes
    clean_branch_t rec(event.Reconstructed);
    clean_branch_t mc(event.MCTrue);

    // run the physics classes
    for( auto& m : physics ) {
        m->ProcessEvent(event, manager);
    }
}

void PhysicsManager::SaveEvent(slowcontrol::event_t buffered_event, const physics::manager_t& manager)
{
    auto& event = buffered_event.Event;

    if(manager.saveEvent || buffered_event.Save || event->SavedForSlowControls) {
        if(!buffered_event.Save && treeEvents->GetCurrentFile() == nullptr)
            LOG_N_TIMES(1, WARNING) << "Writing treeEvents to memory. Might be a lot of data!";


        // always keep read hits if saving for slowcontrol
        if(!manager.keepReadHits && !buffered_event.Save)
            event->Reconstructed->DetectorReadHits.resize(0);

        // indicate that this event was saved for slowcontrol purposes only
        if(!manager.saveEvent && buffered_event.Save) {
            event->SavedForSlowControls = true;
        }

        treeEventPtr = event.get();
        treeEvents->Fill();
    }
}
