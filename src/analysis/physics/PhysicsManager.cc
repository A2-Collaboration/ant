#include "PhysicsManager.h"

#include "utils/ParticleID.h"
#include "input/DataReader.h"

#include "tree/TSlowControl.h"
#include "tree/TAntHeader.h"
#include "base/Logger.h"

#include "input/detail/SlowcontrolCreator.h"

#include "base/ProgressCounter.h"

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

bool PhysicsManager::InitReaders(PhysicsManager::readers_t readers_)
{
    readers = move(readers_);

    // figure out what set of readers we have
    // we expect maximum one source and several amends
    auto it_reader = readers.begin();
    while(it_reader != readers.end()) {
        if((*it_reader)->IsSource()) {
            if(source != nullptr) {
                LOG(ERROR) << "Found more than one source for events, stop.";
                return false;
            }
            source = move(*it_reader);
            it_reader = readers.erase(it_reader);
        }
        else {
            ++it_reader;
        }
    }
    return true;
}

bool PhysicsManager::TryReadEvent(unique_ptr<data::Event>& event)
{
    bool read_event = false;

    while(!read_event) {
        if(source) {
            bool event_ok = source->ReadNextEvent(*event);
            auto sc = source->ReadNextSlowControl();
            if(!event_ok && !sc)
                return false;
            if(event_ok)
                read_event = true;
            if(sc)
               slowcontrol_mgr.ProcessSlowcontrol(move(sc));
        }

        auto it_reader = readers.begin();
        while(it_reader != readers.end()) {
            bool event_ok = (*it_reader)->ReadNextEvent(*event);
            auto sc = (*it_reader)->ReadNextSlowControl();

            if(!event_ok && !sc) {
                it_reader = readers.erase(it_reader);
                continue;
            }

            if(event_ok)
                read_event = true;
            if(sc)
                slowcontrol_mgr.ProcessSlowcontrol(move(sc));
            ++it_reader;
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
        auto& eventid = event->Reconstructed.Trigger.EventID;

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
    if(!InitReaders(move(readers_)))
        return;

    if(physics.empty()) {
        throw Exception("No Analysis Instances activated. Will not analyse anything.");
    }



    auto slkeys = RequestedKeys(slowcontrol_data);

    VLOG(7) << "Requested Slowcontrol keys";
    for(const auto& key : slkeys) {
        VLOG(7) << key;
    }

    slowcontrol_mgr.SetRequiredKeys(slkeys);



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
            auto event = std_ext::make_unique<data::Event>();
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
              << nEventsProcessed/elapsed_seconds.count() << " Items/s";
}





void PhysicsManager::ProcessEvent(unique_ptr<data::Event> event)
{
    if(particleID) {
        // run particle ID for Reconstructed candidates
        auto& reconstructed = event->Reconstructed;
        for(const auto& cand : reconstructed.Candidates) {
            auto particle = particleID->Process(cand);
            if(particle)
                reconstructed.Particles.AddParticle(particle);
        }
    }

    for( auto& m : physics ) {
        m->ProcessEvent(*event);
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
