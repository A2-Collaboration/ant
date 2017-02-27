#pragma once

#ifndef __CINT__

#include <string>
#include <vector>


#include "mc/database/Query.h"

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


    //-- Options ---
    std::string _outfileName;
    std::vector<double> _energies;
    bool _saveUnstable;
    bool _doBulk;
    TF1 _energyFunction;
    const data::Query::ChannelSelector_t ChannelSelector;

    //-- Output ---
    /// \todo use WrapTFile here
    TFile* _outfile;
    TTree* _data;

    //-- data:
    std::vector<BinContent> _energyBins;

    //-- Tools ---
    TRandom3* _rndEngine;

    void initLUT();
    PReaction* makeReaction(const double energy, const ParticleTypeTreeDatabase::Channel& channel) const;

    /**
     * @brief getRandomReaction
     * @return pointer to randomly picked Pluto reaction from database
     */
    PReaction* getRandomReaction() const;

    /// called in destructor
    virtual void finish() const;

public:

    Cocktail(const std::string& outfile,
             const std::vector<double>& energies,
             bool saveUnstable = 0, bool doBulk = 1,
             const std::string& energyDistribution = "1.0 / x",
             const data::Query::ChannelSelector_t& selector
                        = data::Query::GetSelector(data::Query::Selection::All));

    virtual unsigned long Sample(const unsigned long &nevts) const override;

    virtual ~Cocktail(){ finish(); }

};





} //mc
} //simulation
} //ant

#endif // __CINT__

