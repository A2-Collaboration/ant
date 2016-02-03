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

PhysicsManager::PhysicsManager(volatile bool* running_) :
    physics(),
    running(running_)
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

void PhysicsManager::InitReaders(PhysicsManager::readers_t readers_)
{
    readers = move(readers_);

    // figure out what set of readers we have
    // we expect maximum one source and several amends
    auto it_reader = readers.begin();
    while(it_reader != readers.end()) {
        if((*it_reader)->IsSource()) {
            if(source != nullptr)
                throw Exception("Found more than one source in given readers");
            source = move(*it_reader);
            it_reader = readers.erase(it_reader);
        }
        else {
            ++it_reader;
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

    bool finished_reading = false;

    ProgressCounter progress;
    while(true) {
        if(finished_reading)
            break;

        if(!running) {
            VLOG(3) << "Stop requested";
            finished_reading = true;
            break;
        }

        // read new event

        auto event = std_ext::make_unique<TEvent>();
        if(!TryReadEvent(event)) {
            VLOG(5) << "No more events to read, finish.";
            finished_reading = true;
            break;
        }

        slowcontrol_mgr->ProcessEvent(move(event));

        while(auto buffered_event = slowcontrol_mgr->PopEvent()) {

            if(nEventsProcessed == maxevents) {
                VLOG(3) << "Reached max Events (" << maxevents;
                finished_reading = false;
                break;
            }

            if(!running)
                break;

            const auto& eventid = buffered_event->Reconstructed->ID;
            if(nEventsProcessed==0)
                firstID = eventid;
            lastID = eventid;

            logger::DebugInfo::nProcessedEvents = nEventsProcessed;
            ProcessEvent(move(buffered_event));

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
    LOG(INFO) << "Processed " << nEventsProcessed << " events, speed "
              << nEventsProcessed/elapsed_seconds.count() << " event/s";

    const auto nEventsSaved = treeEvents->GetEntries();
    if(nEventsSaved == 0)
        delete treeEvents;
    else if(treeEvents->GetCurrentFile() != nullptr) {
        treeEvents->Write();
        LOG(INFO) << "Wrote " << nEventsSaved << " treeEvents: "
                  << (double)treeEvents->GetTotBytes()/(1 << 20) << " MB (uncompressed), "
                  << (double)treeEvents->GetTotBytes()/nEventsSaved << " bytes/event";
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


    auto it_reader = readers.begin();
    while(it_reader != readers.end()) {

        if((*it_reader)->ReadNextEvent(*event)) {
            ++it_reader;
            event_read = true;
        }
        else {
            it_reader = readers.erase(it_reader);
        }
    }

    return event_read;
}

void PhysicsManager::ProcessEvent(std::unique_ptr<TEvent> event)
{


    if(particleID && event->Reconstructed) {
        // run particle ID for Reconstructed candidates
        // but only if there are no identified particles present yet
        /// \todo implement flag to force particle ID again?
        TEventData& recon = *event->Reconstructed;
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

    bool clean_reconstructed = false;
    if(!event->Reconstructed) {
        event->Reconstructed = std_ext::make_unique<TEventData>();
        clean_reconstructed = true;
    }
    bool clean_mctrue = false;
    if(!event->MCTrue) {
        event->MCTrue = std_ext::make_unique<TEventData>();
        clean_mctrue = true;
    }

    // run the physics classes
    processmanager.Reset();
    for( auto& m : physics ) {
        m->ProcessEvent(*event, processmanager);
    }

    // remove the temporary empty branches again
    if(clean_reconstructed)
        event->Reconstructed = nullptr;
    if(clean_mctrue)
        event->MCTrue = nullptr;

    if(processmanager.saveEvent) {
        if(treeEvents->GetCurrentFile() == nullptr)
            LOG_N_TIMES(1, WARNING) << "Writing treeEvents to memory. Might be a lot of data!";
        if(!processmanager.keepReadHits)
            event->Reconstructed->DetectorReadHits.resize(0);
        treeEventPtr = event.get();
        treeEvents->Fill();
    }

}
