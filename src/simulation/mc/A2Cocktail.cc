#include "A2Cocktail.h"

#include <iostream>

#include "PlutoExtensions.h"

#include "TH1D.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wwrite-strings"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wvla-extension"
#endif

#include "PDecayManager.h"
#include "PDecayChannel.h"
#include "PReaction.h"
#include "PParticle.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif
#pragma GCC diagnostic pop

using namespace std;
using namespace ant::simulation::mc;


void FixedEnergyCocktail::_makeDecays()
{
    A2ChannelManager a2man;
    _primary_decays = a2man.GenerateDecays(_E);
}

FixedEnergyCocktail::FixedEnergyCocktail(const string& outfile, double Emin, double Emax, bool bulk, bool stable):
    _E(( Emin + Emax ) / 2.0),
    _stable(stable),
    ofile(outfile)
{

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

unsigned long FixedEnergyCocktail::Sample(const unsigned long& nevts) const
{
    return _pdm->Loop(nevts,0,strdup(ofile.c_str()) ,_stable,0,1,0,1);
    //                                             ...,stable,obs,vertex,ascii,random)
    //                                                                   files  order
}





PReaction *A2Cocktail::GetRandomReaction() const
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

            if ( xsection != 0)
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

PReaction *A2Cocktail::makeReaction(const double& energy, const string &particles) const
{
    // convert database entry t pluto-decay-string: g + p --> p + product1 + ...
    string reactionstring = "p " + particles;

    PReaction* reaction = new PReaction(energy,                         // beam momentum = photon engry
                                        strdup("g"),strdup("p"),        // beam,target
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
    UpdatePluteDataBase();
}




unsigned long A2Cocktail::Sample(const unsigned long &nevts) const
{
    unsigned long errors(0);

    for ( unsigned long evt = 0 ; evt < nevts ; ++evt){
        errors += 1 - GetRandomReaction()->Loop(1,0,0); // reminder: Loop(numEvents,weightFlag,verbose)....
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
