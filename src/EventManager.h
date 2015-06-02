#ifndef EVENTMANAGER_H
#define EVENTMANAGER_H

// GoAT headers
#include "GTreeManager.h"
#include "GTree.h"
#include "GTreeTrack.h"
#include "GTreeParticle.h"
#include "GTreeMeson.h"
#include "GTreeDetectorHits.h"
#include "GTreeEventParameters.h"
#include "GTreeA2Geant.h"
#include "GTreePluto.h"

// Ant
#include "ParticleType.h"
#include "Particle.h"
#include "Event.h"
#include "AntPhysics.h"

// std
#include <stdexcept>
#include <memory>
#include <stdexcept>


class PStaticData;

namespace ant {

/**
 * @brief data_check_exception is thrown if a data input quality check fails
 */
class data_check_exception : public std::exception {
protected:
    std::string msg;

public:
    data_check_exception(const std::string& message): msg(message) {}
    const char *what() const throw() { return msg.c_str(); }
};


/**
 * @brief The EventManager class
 */
class EventManager: public GTreeManager {
public:
    using PhysicsList = std::list<ant::Physics*>;

protected:

    unsigned int maxevents;

    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();

    PhysicsList physics;

    void RunPhysics(const ant::Event& event);


    void CopyParticles(GTreeParticle* tree, const ant::ParticleTypeDatabase::Type& type, ant::Event& event);
    void CopyTracks(GTreeTrack* tree, Event& event);
    void CopyPlutoParticles(GTreePluto* tree, Event &event);
    void CopyTaggerHits(Event& event);
    void CopyTriggerInfo(GTreeTrigger* tree, Event& event);

    void checkMCIDs();

    PStaticData* pluto_database;

public:
    EventManager();
    virtual ~EventManager();
    virtual Bool_t  Init(const char* configfile);

    void AddPhysics(ant::Physics* phys) { physics.push_back(phys); }
    void SetMaxEvents(unsigned int max) { maxevents = max; }

    void Finish();

};
}

#endif
