#include "PlutoGenerator.h"


#include "PlutoExtensions.h"

#include "PReaction.h"
#include "PParticle.h"

using namespace std;
using namespace ant::mc;
using namespace ant::mc::pluto;

PlutoGenerator::PlutoGenerator(const bool updateDataBase)
{
    if (updateDataBase)
        UpdatePlutoDataBase();
}


Cocktail::Cocktail(const string& outfile,
                   const std::vector<double>& energies,
                   bool saveUnstable, bool doBulk,
                   const int verbosity,
                   const string& energyDistribution,
                   const data::Query::ChannelSelector_t& selector):
    _fileOutput(outfile),
    _energies(energies),
    _settings(saveUnstable,doBulk),
    ChannelSelector(selector)
{
    sort(_energies.begin(), _energies.end());
    _energyFunction = TF1("beamEnergy",energyDistribution.c_str(),_energies.front(),_energies.back());
#ifdef PLUTO_GLOBAL_H
    pluto_global::verbosity = verbosity;
#endif
    initLUT();
}

void Cocktail::initLUT()
{
    // helpers:
    double acc_E = 0;

    // -- Init outputfile and Tree --
    _data = _fileOutput.CreateInside<TTree>("data","Event data");

    // -- Init root - random engine ---
    _rndEngine = new TRandom3(0);

    for(double energy : _energies)
    {
        BinContent currentBin;
        // -- Energy in MeV --
        currentBin.Energy = energy;

        // -- Statistics for Energy --
        //    p(E) = f(E) * totalXsection(E)
        acc_E += _energyFunction(currentBin.Energy) * data::Query::TotalXsection(currentBin.Energy,ChannelSelector);
        currentBin.AccProbability = acc_E;

        // -- Statistics for Channels in this energy bin --
        //    Generate accumulated probabilities with PReaction(tm)
        //    for all channels in at current energy
        double acc_prob_channels = 0;
        for ( auto& product: data::Query::GetProductionChannels(ChannelSelector))
        {
            double xsection = data::Query::Xsection(product, currentBin.Energy);

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



PReaction *Cocktail::makeReaction(const double energy,
                                  const ParticleTypeTreeDatabase::Channel& channel) const
{
    auto reaction = PlutoFactory::MakeFixedEnergyReaction(
                energy,
                channel,
                _data,
                _settings);

    return reaction;
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
