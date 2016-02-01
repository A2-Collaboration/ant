#include "PhysicsManager.h"

#include "utils/ParticleID.h"
#include "input/DataReader.h"

#include "tree/TSlowControl.h"
#include "tree/TAntHeader.h"
#include "base/Logger.h"

#include "input/slowcontrol/SlowControlCreator.h"

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

bool PhysicsManager::TryReadEvent(TEventPtr& event)
{
    bool read_event = false;

    while(!read_event) {
        if(source) {
            if(source->ReadNextEvent(*event)) {
                read_event = true;
                slowcontrol_mgr.ProcessSlowControls(*event);
            }
            else {
                return false;
            }
        }

        auto it_reader = readers.begin();
        while(it_reader != readers.end()) {

            if((*it_reader)->ReadNextEvent(*event)) {
                read_event = true;
                slowcontrol_mgr.ProcessSlowControls(*event);
                ++it_reader;
            }
            else {
                it_reader = readers.erase(it_reader);
            }
        }

        if(!source && readers.empty())
            return false;
    }
    return true;
}

void PhysicsManager::ProcessEventBuffer(long long maxevents)
{
    // if running is already false,
    // flush the buffer no matter what...
    bool flush = !running;
    if(flush)
        VLOG(5) << "Flushing " << eventbuffer.size() << " events from eventbuffer";

    TID runUntil = slowcontrol_mgr.UpdateSlowcontrolData(slowcontrol_data);

    if(slowcontrol_mgr.hasRequests() && runUntil.IsInvalid())
        return;

    while(!eventbuffer.empty()) {
        // running might change to false here in this loop
        if(!running && !flush)
            return;

        if(nEventsProcessed>=maxevents)
            return;

        auto& event = eventbuffer.front();
        auto& eventid = event->Reconstructed->ID;

        if(slowcontrol_mgr.hasRequests() && (eventid > runUntil))
            break;

        ProcessEvent(move(event));
        eventbuffer.pop();
        if(nEventsProcessed==0)
            firstID = eventid;
        lastID = eventid;
        nEventsProcessed++;
    }
}



void PhysicsManager::ReadFrom(
        list<unique_ptr<input::DataReader> > readers_,
        long long maxevents)
{
    InitReaders(move(readers_));

    if(physics.empty())
        throw Exception("No analysis instances activated. Cannot not analyse anything.");

    // prepare slowcontrol
    auto slkeys = RequestedKeys(slowcontrol_data);
    VLOG(7) << "Requested Slowcontrol keys";
    for(const auto& key : slkeys) {
        VLOG(7) << key;
    }
    slowcontrol_mgr.SetRequiredKeys(slkeys);


    // prepare output of TEvents
    treeEvents = new TTree("treeEvents","TEvent data");
    treeEventPtr = nullptr;
    treeEvents->Branch("data", addressof(treeEventPtr));


    chrono::time_point<std::chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    bool finished_reading = false;
    long long nEventsRead = 0;

    ProgressCounter progress;
    while(true) {
        if(finished_reading)
            break;

        do
        {
            if(!running || nEventsRead>=maxevents) {
                VLOG(5) << "End of reading requested";
                finished_reading = true;
                break;
            }
            auto event = std_ext::make_unique<TEvent>();
            if(!TryReadEvent(event)) {
                VLOG(5) << "No more events to read, finish.";
                finished_reading = true;
                break;
            }
            if(eventbuffer.size()>20000) {
                throw Exception("SlowControl buffering reached maximum size without becoming complete.");
            }

            eventbuffer.emplace(move(event));
            nEventsRead++;

            if(progressUpdates && source && progress.Update(source->PercentDone())) {
                LOG(INFO) << progress;
            }

        }
        while(!slowcontrol_mgr.isComplete());

        if(slowcontrol_mgr.isComplete()) {
            VLOG(5) << "Slowcontrol set complete. Processing event buffer.";
            ProcessEventBuffer(maxevents);
        } else {
            VLOG(5) << "Finished reading but slowcontrol block not complete. Dropping remaining " << eventbuffer.size() << " events.";
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

void PhysicsManager::ProcessEvent(std::unique_ptr<TEvent> event)
{
    if(particleID && event->Reconstructed) {
        // run particle ID for Reconstructed candidates
        // but only if there are no identified particles present yet
        /// \todo implement flag to force particle ID again?
        TEvent::Data& recon = *event->Reconstructed;
        if(recon.Particles.GetAll().empty()) {
            for(const auto& cand : recon.Candidates) {
                auto particle = particleID->Process(cand);
                if(particle)
                    recon.Particles.Add(particle);
            }
        }
    }

    // ensure that physics classes always
    // have at least empty TEvent::Data branches MCTrue and Reconstructed

    bool clean_reconstructed = false;
    if(!event->Reconstructed) {
        event->Reconstructed = std_ext::make_unique<TEvent::Data>();
        clean_reconstructed = true;
    }
    bool clean_mctrue = false;
    if(!event->MCTrue) {
        event->MCTrue = std_ext::make_unique<TEvent::Data>();
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

void PhysicsManager::ShowResults()
{
    for(auto& p : physics) {
        p->ShowResult();
    }
}

void PhysicsManager::EnableProgressUpdates(bool updates)
{
    progressUpdates = updates;
}
