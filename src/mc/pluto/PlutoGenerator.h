#pragma once

#ifndef __CINT__

#include <string>
#include <vector>

#include "base/WrapTFile.h"

#include "mc/database/Query.h"
#include "PlutoFactory.h"



// ROOOOOT
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TRandom3.h"

class PReaction;
class PParticle;
class PPlutoBulkDecay;

namespace ant
{
namespace mc
{
namespace pluto
{


class PlutoGenerator{
public:
    /**
     * @brief Sample: Main loop
     * @param nevts
     * @return number of failed events
     */
    virtual unsigned long Sample(const unsigned long& nevts) const=0;

protected:

    /// Ensure every derived class loads our plotu extensions as default!
    PlutoGenerator( const bool updateDataBase = true);

    ~PlutoGenerator() = default;
};


/**
 * @brief The A2Cocktail
 */
class Cocktail: public PlutoGenerator
{
private:
    /**
     * @brief The BinContent class is a container for all necessary information in each energy channel.
     *        The data is provided by A2ChannelManager
     */
    struct BinContent{
        struct Reaction_t{
            const double AccProb;
            PReaction*   PlutoReaction;
            Reaction_t(double accProb, PReaction* pReaction):
                AccProb(accProb),
                PlutoReaction(pReaction){}
        };

        double Energy; // in MeV
        double AccProbability;
        std::vector<std::string> DecayProducts;
        std::vector<Reaction_t>  ReactionLUT;
    };


    WrapTFileOutput _fileOutput;
    //-- Options ---
    std::vector<double> _energies;
    ReactionSettings_t _settings;
    TF1 _energyFunction;
    const data::Query::ChannelSelector_t ChannelSelector;


    TTree* _data;

    //-- data:
    std::vector<BinContent> _energyBins;

    //-- Tools ---
    TRandom3* _rndEngine;

    void initLUT();
    PReaction* makeReaction(const double energy,
                            const ParticleTypeTreeDatabase::Channel& channel) const;

    /**
     * @brief getRandomReaction
     * @return pointer to randomly picked Pluto reaction from database
     */
    PReaction* getRandomReaction() const;

public:

    Cocktail(const std::string& outfile,
             const std::vector<double>& energies,
             bool saveUnstable = 0, bool doBulk = 1,
             const int verbosity = 0,
             const std::string& energyDistribution = "1.0 / x",
             const data::Query::ChannelSelector_t& selector
                        = data::Query::GetSelector(data::Query::Selection::gpBeamTarget));

    virtual unsigned long Sample(const unsigned long &nevts) const override;

    virtual ~Cocktail(){}

};





} //mc
} //simulation
} //ant

#endif // __CINT__

