#pragma once

#ifndef __CINT__

#include <string>
#include <vector>


#include "A2Channels.h"

// ROOOOOT
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TRandom3.h"

class PDecayManager;
class PDecayChannel;
class PReaction;
class PParticle;
class PPlutoBulkDecay;

namespace ant
{
namespace simulation
{
namespace mc
{


class ManagedPlutoReaction{
public:
    /**
     * @brief Sample: Main loop
     * @param nevts
     * @return number of failed events
     */
    virtual unsigned long Sample(const unsigned long& nevts) const=0;
};


/**
 * @brief The A2Cocktail
 */
class A2Cocktail: public ManagedPlutoReaction
{
private:
    /**
     * @brief The BinContent class is a container for all necessary information in each energy channel.
     *        The data is provided by A2ChannelManager
     */
    class BinContent{
    public:
        double Energy; // in GeV
        double AccProbability;
        std::vector<std::string> DecayProducts;
        std::vector<std::pair<double,PReaction*>> Channellist;
    };

    //-- Options ---
    std::string _outfileName;
    std::vector<double> _energies;
    bool _saveUnstable;
    bool _doBulk;
    TF1 _energyFunction;

    //-- Output ---
    TFile* _outfile;
    TTree* _data;

    //-- data:
    std::vector<BinContent> _energyBins;

    //-- Tools ---
    TRandom3* _rndEngine;

    void init(std::vector<std::string> filelist);
    PReaction* makeReaction(const double &energy, const std::string& particles) const;

public:

    /**
     * @brief constuctor note: Provide Enegies in GeV
     *
     * @return pointer to randomly picked Pluto reaction from database
     */
    A2Cocktail(const std::string& outfile,
               const std::vector<double>& energies,
               bool saveUnstable = 0, bool doBulk = 1,
               std::vector<std::string> filenames = {},
               const std::string& energyDistribution = "1.0 / x" );

    /**
     * @brief getRandomReaction
     * @return pointer to randomly picked Pluto reaction from database
     */
    PReaction* GetRandomReaction() const;

    virtual unsigned long Sample(const unsigned long &nevts) const override;
    virtual void Finish() const;
    //virtual ~A2Cocktail(){Finish();} maybe better, not sure yet
};


/**
 * @brief The A2OldCocktail class generates MC events for a fixed energy, using the Pluto class PDecayManager
 */
// @ Andi bitte nicht loeschen!!!!!
class FixedEnergyCocktail: public ManagedPlutoReaction
{
private:

    double _E;
    bool _stable;
    std::string ofile;

    PDecayManager* _pdm;
    PParticle* _fusion;
    PDecayChannel* _primary_decays;
    PPlutoBulkDecay* _bulkdecay;

    void _makeDecays();

public:
    FixedEnergyCocktail(const std::string& outfile, double Emin, double Emax, bool bulk, bool stable);

    virtual unsigned long Sample(const unsigned long& nevts) const override;

};


} //mc
} //simulation
} //ant

#endif // __CINT__

