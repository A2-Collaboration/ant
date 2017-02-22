#include "A2Cocktail.h"

#include <iostream>

#include "PlutoExtensions.h"

#include "TH1D.h"

#include "PDecayManager.h"
#include "PDecayChannel.h"
#include "PReaction.h"
#include "PParticle.h"

using namespace std;
using namespace ant::mc::pluto;




PReaction* A2Cocktail::getRandomReaction() const
{
    double rndEnergyBinValue = _rndEngine->Rndm() * _energyBins.back().AccProbability;

    for (auto& eBin: _energyBins)
    {
        if ( rndEnergyBinValue <= eBin.AccProbability )
        {
            double rndChannelValue = _rndEngine->Rndm() * eBin.Channellist.back().first;
            for (auto& channel: eBin.Channellist)
            {
                if ( rndChannelValue <= channel.first ){
                    return channel.second;
                }
            }
        }
    }
    return nullptr;
}

void A2Cocktail::init(vector<string> filelist)
{
    // helpers:
    double acc_E = 0;
    BinContent currentBin;
    A2ChannelManager a2man(filelist);

    // -- Init outputfile and Tree --
    _outfile = new TFile(string(_outfileName + ".root").c_str(),"recreate");
    _data = new TTree("data","Event data");

    // -- Init root - random engine ---
    _rndEngine = new TRandom3(0);

    for(double energy : _energies)
    {
        // -- Channel product names --
        currentBin.DecayProducts = a2man.GetChannels();

        // -- Energy in GeV --
        currentBin.Energy = energy;

        // -- Statistics for Energy --
        //    p(E) = f(E) * totalXsection(E)
        acc_E += _energyFunction(currentBin.Energy) * a2man.TotalXsection(currentBin.Energy);
        currentBin.AccProbability = acc_E;

        // -- Statistics for Channels in this energy bin --
        //    Generate accumulated probabilities with PReaction(tm)
        //    for all channels in at current energy
        double acc_prob_channels = 0;
        for ( auto& product: currentBin.DecayProducts)
        {
            double xsection = a2man.Xsection(product, currentBin.Energy);

            if ( xsection > 0)
            {
                acc_prob_channels += xsection;
                pair<double,PReaction*> accprobReactionPair;

                accprobReactionPair.first  = acc_prob_channels;
                accprobReactionPair.second = makeReaction(currentBin.Energy,product);
                currentBin.Channellist.push_back(accprobReactionPair);
            }
        }
        // -- fill --
        _energyBins.push_back(currentBin);
        currentBin.Channellist.clear();
    }
}

/// TODO allow different target particles here
PReaction *A2Cocktail::makeReaction(const double& energy, const string &outGoingParticles, const string& targetParticle) const
{
    // convert database entry t pluto-decay-string: g + p --> p + product1 + ...
    //assume only reactions g p -> X p
    /// TODO allow targetParticle not in outgoing particles (EG g p -> sigma+ k0)
    string reactionstring = targetParticle + " " + outGoingParticles;

    PReaction* reaction = new PReaction(energy,                         // beam momentum = photon engry
                                        strdup("g"),strdup(targetParticle.c_str()),        // beam,target
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
    return reaction;
}

A2Cocktail::A2Cocktail(const string& outfile,
                       const std::vector<double>& energies,
                       bool saveUnstable, bool doBulk,
                       std::vector<string> filenames,
                       const string& energyDistribution):
    _outfileName(outfile),
    _energies(energies),
    _saveUnstable(saveUnstable),
    _doBulk(doBulk)
{
    sort(_energies.begin(), _energies.end());
    _energyFunction = TF1("beamEnergy",energyDistribution.c_str(),_energies.front(),_energies.back());
    init(filenames);
    UpdatePlutoDataBase();
}




unsigned long A2Cocktail::Sample(const unsigned long &nevts) const
{
    unsigned long errors(0);

    for ( unsigned long evt = 0 ; evt < nevts ; ++evt){
        errors += 1 - getRandomReaction()->Loop(1,0,0); // reminder: Loop(numEvents,weightFlag,verbose)....
    }
    return errors;
}

void A2Cocktail::Finish() const
{
    if (_outfile)
    {
        _outfile->Write();
        _outfile->Close();
    }
}
