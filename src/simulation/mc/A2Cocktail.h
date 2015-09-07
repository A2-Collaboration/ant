/**
 * @file A2ocktail.h
 * @author Martin Wolfes
 * @date June 2015
 *
 */

#ifndef _A2COCKTAIL_H
#define _A2COCKTAIL_H

#ifndef __CINT__



#include <string>
#include <vector>

// pluto++
#pragma GCC diagnostic ignored "-Wwrite-strings"
#include "PDecayManager.h"
#include "PDecayChannel.h"
#include "PReaction.h"
#include "PParticle.h"

#include "A2Channels.h"

// ROOOOOT
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TRandom3.h"


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
class A2Cocktail: public ManagedPlutoReaction{
private:
    /**
     * @brief The BinContent class is a container for all necessary information in each energy channel.
     *        The data is provided by A2ChannelManager
     */
    class BinContent{
    public:
        double Energy;
        double AccProbability;
        std::vector<std::string> DecayProducts;
        std::vector<std::pair<double,PReaction*>> Channellist;
    };

    //-- Options ---
    std::string _outfileName;
    double _Emin;
    double _Emax;
    double _numEnergyBins;
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
    A2Cocktail(const std::string& outfile,
               const double& Emin, const double& Emax,
               const unsigned int numEnergyBins,
               bool saveUnstable = 0, const bool doBulk = 1,
               std::vector<std::string> filenames = {},
               const std::string& energyDistribution = "1.0 / x" ):
        _outfileName(outfile),
        _Emin(Emin), _Emax(Emax),
        _numEnergyBins(numEnergyBins),
        _saveUnstable(saveUnstable),
        _doBulk(doBulk),
        _energyFunction("beamEnergy",strdup(energyDistribution.c_str()),Emin,Emax)
    {
        init(filenames);
    }

    /**
     * @brief getRandomReaction
     * @return pointer to randomly picked Pluto reaction from database
     */
    PReaction* GetRandomReaction() const;

    virtual unsigned long Sample(const unsigned long &nevts) const;
    virtual void Finish() const;
    //virtual ~A2Cocktail(){Finish();} maybe better, not sure yet
};


/**
 * @brief The A2OldCocktail class generates MC events for a fixed energy, using the Pluto class PDecayManager
 */
class A2OldCocktail{
private:

    double _E;
    bool _stable;

    PDecayManager* _pdm;
    PParticle* _fusion;
    PDecayChannel* _primary_decays;
    PPlutoBulkDecay* _bulkdecay;

    void _makeDecays();

public:
    A2OldCocktail(double Emin, double Emax, bool bulk, bool stable):
        _E(( Emin + Emax ) / 2.0), _stable(stable){

        _fusion = new PParticle( PParticle(1,_E) + PParticle(14));
        _pdm = new PDecayManager();

        _makeDecays();

        _pdm->InitReaction(_fusion,_primary_decays);

        if(bulk){
            _bulkdecay = new PPlutoBulkDecay();
            _bulkdecay->SetRecursiveMode(1);
            _bulkdecay->SetTauMax(0.001);
            _pdm->AddBulk(_bulkdecay);
        }
    }
    int Sample(unsigned int n_evts, const std::string& ofile) {
        return _pdm->Loop(n_evts,0,strdup(ofile.c_str()) ,_stable,0,1,0,1);
        //                                             ...,stable,obs,vertex,ascii,random)
        //                                                                   files  order
    }

};

#endif 

#endif 
