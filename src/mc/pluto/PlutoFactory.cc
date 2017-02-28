#include "PlutoFactory.h"

#include "mc/database/Query.h"
#include "PReaction.h"

#include "base/Logger.h"

#include "TTree.h"

using namespace std;
using namespace ant::mc;
using namespace ant::mc::pluto;


PPlutoBulkDecay* addRecursiveBulkDecay(PReaction* reaction, const double tauMax)
{
    PPlutoBulkDecay* bulkdecay = new PPlutoBulkDecay();
    bulkdecay->SetRecursiveMode(1);
    bulkdecay->SetTauMax(tauMax);
    reaction->AddBulk(bulkdecay);
    return bulkdecay;
}

void preHeat(PReaction* reaction, const int amount)
{
    reaction->Preheating(amount);
}

PReaction*PlutoFactory::MakeFixedEnergyReaction(const double beamMomentum,
                                                const ant::ParticleTypeTreeDatabase::Channel& channel,
                                                TTree* outTree,
                                                const ReactionSettings_t& settings)
{
    const auto reactionstring = data::Query::GetPlutoProductString(channel);

    const auto beamstring     = data::Query::GetPlutoBeamString(channel);
    const auto targetstring   = data::Query::GetPlutoTargetString(channel);

    if (outTree->GetCurrentFile() == nullptr)
    {
        throw runtime_error("No file attached to pluto-tree.");
    }

    PReaction* reaction = new PReaction(beamMomentum * MeVtoGeV, // beam momentum (pluto uses in GeV)
                                        strdup(beamstring.c_str()),
                                        strdup(targetstring.c_str()),
                                        strdup(reactionstring.c_str()),
                                        strdup(""),   // output - filename left empty to prevent pluto from filehandling
                                        settings.SaveUnstable,
                                        settings.UnusedPReactionFlag,
                                        settings.CalculateVertices,
                                        settings.AsciiOutPut,
                                        outTree);

    if (settings.DoRecusiveBulkDecay)
    {
        addRecursiveBulkDecay(reaction, settings.MaxTauForBulkDecay);
    }

    preHeat(reaction,settings.preheatedEvents);

    return reaction;
}
