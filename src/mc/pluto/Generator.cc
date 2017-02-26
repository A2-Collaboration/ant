#include "Generator.h"

#include <iostream>

#include "PlutoExtensions.h"

#include "TH1D.h"

#include "PDecayManager.h"
#include "PDecayChannel.h"
#include "PReaction.h"
#include "PParticle.h"

using namespace std;
using namespace ant::mc::pluto;




PReaction* Cocktail::getRandomReaction() const
{
    double rndEnergyBinValue = _rndEngine->Rndm() * _energyBins.back().AccProbability;

    for (auto& eBin: _energyBins)
    {
        if ( rndEnergyBinValue <= eBin.AccProbability )
        {
            double rndChannelValue = _rndEngine->Rndm() * eBin.ReactionLUT.back().AccProb;
            for (auto& reactrion: eBin.ReactionLUT)
            {
                if ( rndChannelValue <= reactrion.AccProb ){
                    return reactrion.PlutoReaction;
                }
            }
        }
    }
    return nullptr;
}

void Cocktail::init()
{
    // helpers:
    double acc_E = 0;
    A2ChannelManager chMan;

    // -- Init outputfile and Tree --
    _outfile = new TFile(string(_outfileName + ".root").c_str(),"recreate");
    _data = new TTree("data","Event data");

    // -- Init root - random engine ---
    _rndEngine = new TRandom3(0);

    for(double energy : _energies)
    {           
        BinContent currentBin;
        // -- Energy in GeV --
        currentBin.Energy = energy;

        // -- Statistics for Energy --
        //    p(E) = f(E) * totalXsection(E)
        acc_E += _energyFunction(currentBin.Energy) * chMan.TotalXsection(currentBin.Energy);
        currentBin.AccProbability = acc_E;

        // -- Statistics for Channels in this energy bin --
        //    Generate accumulated probabilities with PReaction(tm)
        //    for all channels in at current energy
        double acc_prob_channels = 0;
        for ( auto& product: chMan.GetChannels())
        {
            double xsection = chMan.Xsection(product, currentBin.Energy);

            if ( xsection > 0)  // make sure channel is available
            {
                acc_prob_channels += xsection;

                currentBin.ReactionLUT.emplace_back(BinContent::Reaction_t(acc_prob_channels,
                                                                           makeReaction(currentBin.Energy,product)));
            }
        }
        // -- fill --
        _energyBins.emplace_back(currentBin);
    }
}

/// TODO allow different target particles here
PReaction *Cocktail::makeReaction(const double energy, const ParticleTypeTreeDatabase::Channel& channel) const
{
    // convert database entry t pluto-decay-string: g + p --> p + product1 + ...
    //assume only reactions g p -> X p
    /// TODO allow targetParticle not in outgoing particles (EG g p -> sigma+ k0)
    ///
    A2ChannelManager chMan;

    auto reactionstring = chMan.GetPlutoProductString(channel);

    auto beamstring     = chMan.GetBeam(channel);
    auto targetstring   = chMan.GetTarget(channel);


    PReaction* reaction = new PReaction(energy,                         // beam momentum = photon engry
                                        strdup(beamstring.c_str()),strdup(targetstring.c_str()),        // beam,target
                                        strdup(reactionstring.c_str()),
                                        strdup(_outfileName.c_str()),   // output - filename
                                        _saveUnstable,0,1,0,            // pluto - flags
                                        _data);                         // pointer to output tree

    // -- Init bulk-interface if needed ---
    if (_doBulk)
    {
        PPlutoBulkDecay* bulkdecay = new PPlutoBulkDecay();
        bulkdecay->SetRecursiveMode(1);
        bulkdecay->SetTauMax(0.001);
        reaction->AddBulk(bulkdecay);
    }
    reaction->Print();
    return reaction;
}

Cocktail::Cocktail(const string& outfile,
                       const std::vector<double>& energies,
                       bool saveUnstable, bool doBulk,
                       std::vector<string>,
                       const string& energyDistribution):
    _outfileName(outfile),
    _energies(energies),
    _saveUnstable(saveUnstable),
    _doBulk(doBulk)
{
    sort(_energies.begin(), _energies.end());
    _energyFunction = TF1("beamEnergy",energyDistribution.c_str(),_energies.front(),_energies.back());
    init();
    UpdatePlutoDataBase();
}




unsigned long Cocktail::Sample(const unsigned long &nevts) const
{
    unsigned long errors(0);

    for ( unsigned long evt = 0 ; evt < nevts ; ++evt){
        errors += 1 - getRandomReaction()->Loop(1,0,0); // reminder: Loop(numEvents,weightFlag,verbose)....
                                                                                        //weight does nothing??!!
    }
    return errors;
}

void Cocktail::Finish() const
{
    if (_outfile)
    {
        _outfile->Write();
        _outfile->Close();
    }
}
