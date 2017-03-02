#pragma once

#ifndef __CINT__

#include "base/ParticleTypeTree.h"
#include "TF1.h"

class PReaction;
class TTree;

namespace ant
{
namespace mc
{
namespace pluto
{

struct ReactionSettings_t{
    bool      SaveUnstable;
    bool      DoRecusiveBulkDecay;
    bool      CalculateVertices;


    // these Options should not be changed
    static constexpr int preheatedEvents          = 1000;
    static constexpr bool     AsciiOutPut         = false;
    static constexpr bool     UnusedPReactionFlag = false;
    static constexpr double   MaxTauForBulkDecay  = 0.001;

    ReactionSettings_t(const bool     saveUnstable         = true,
                       const bool     doRecursiveBulkDecay = true,
                       const bool     calculateVertices    = true
                     ):
        SaveUnstable(saveUnstable),
        DoRecusiveBulkDecay(doRecursiveBulkDecay),
        CalculateVertices(calculateVertices){}
};

/*
struct PlutoBeamSettings_t{
    const std::function<double(const double)> MomentumFunction;
    const std::function<double(const double)> AngularSmearing;
    PlutoBeamSettings_t(const std::function<double(const double)> momentumFunction,
                        const std::function<double(const double)> angularSmearing):
        MomentumFunction(momentumFunction),
        AngularSmearing(angularSmearing){}
};
*/

struct PlutoFactory{

static constexpr double MeVtoGeV = 1.0/1000.0;


static PReaction* MakeFixedEnergyReaction(
        const double beamMomentum,
        const ParticleTypeTreeDatabase::Channel& channel,
        TTree* outTree,
        const ReactionSettings_t& reactionSettings = ReactionSettings_t()
        );
/*
static PReaction* MakeSmearedReaction(
        const PlutoBeamSettings_t&,
        const ParticleTypeTreeDatabase::Channel& channel,
        const ReactionSettings_t& reactionSettings,
        TTree* outTree);
static PReaction* MakeSmearedReaction(
        const PlutoBeamSettings_t&,
        const std::string& plutoBeamName,
        const std::string& plutoTargetName,
        const std::string& plutoDecayString,
        const ReactionSettings_t& reactionSettings,
        TTree* outTree);
        */
};


} //mc
} //simulation
} //ant

#endif // __CINT__

