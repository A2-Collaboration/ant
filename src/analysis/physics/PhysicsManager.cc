#include "PhysicsManager.h"

#include "utils/ParticleID.h"
#include "input/DataReader.h"

#include "tree/TSlowControl.h"
#include "tree/TAntHeader.h"
#include "base/Logger.h"

#include <iomanip>
#include <chrono>

using namespace std;
using namespace ant;
using namespace ant::analysis;

PhysicsManager::PhysicsManager() : physics() {}

PhysicsManager::~PhysicsManager() {}

void PhysicsManager::SetParticleID(std::unique_ptr<utils::ParticleID> pid) {
    particleID = move(pid);
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
        ++it_reader;
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
                slowcontrol[sc->GetKey()].emplace(move(sc));
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
                slowcontrol[sc->GetKey()].emplace(move(sc));
            ++it_reader;
        }

        if(!source && readers.empty())
            return false;
    }
    return true;
}

void PhysicsManager::ProcessEventBuffer(
        long long maxevents,
        bool& running,
        TAntHeader& header
        )
{
    // if running is already false,
    // flush the buffer no matter what...
    bool flush = !running;
    if(flush)
        VLOG(5) << "Flushing " << eventbuffer.size() << " events from eventbuffer";
    while(!eventbuffer.empty()) {
        if(!running && !flush)
            return;
        if(nEventsProcessed>=maxevents)
            return;

        auto& event = eventbuffer.front();
        auto& eventid = event->Reconstructed().TriggerInfos().EventID();
        ProcessEvent(move(event));
        eventbuffer.pop();
        if(nEventsProcessed==0)
            header.FirstID = eventid;
        header.LastID = eventid;
        nEventsProcessed++;
    }
}



void PhysicsManager::ReadFrom(
        std::list<std::unique_ptr<input::DataReader> > readers_,
        long long maxevents,
        bool& running,
        TAntHeader& header
        )
{
    if(!InitReaders(move(readers_)))
        return;

    if(physics.empty()) {
        LOG(WARNING) << "No Analysis Instances activated. Will not analyse anything.";
    }


    chrono::time_point<std::chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    bool finished_reading = false;
    long long nEventsRead = 0;
    while(true) {
        if(finished_reading)
            break;

        while(eventbuffer.size()<20000) {
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
            eventbuffer.emplace(move(event));
            nEventsRead++;
        }

        ProcessEventBuffer(maxevents, running, header);
    }

    VLOG(5) << "First EventId processed: " << header.FirstID;
    VLOG(5) << "Last  EventId processed: " << header.LastID;

    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    LOG(INFO) << "Processed " << nEventsProcessed << " events, speed "
              << nEventsProcessed/elapsed_seconds.count() << " Items/s";
}





void PhysicsManager::ProcessEvent(unique_ptr<data::Event> event)
{
    if(particleID) {
        // run particle ID for Reconstructed candidates
        auto& reconstructed = event->Reconstructed();
        for(const auto& cand : reconstructed.Candidates()) {
            auto particle = particleID->Process(cand);
            if(particle)
                reconstructed.Particles().AddParticle(particle);
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
