#pragma once

#ifndef __CINT__

#include <map>
#include <list>
#include <vector>
#include <string>

#include "ProductionChannel.h"

// ROOT
#include "Math/Interpolator.h"
#include "TH1.h"


class PDecayChannel;

namespace ant
{
namespace mc
{
namespace pluto
{


struct ParticleData
{
    std::vector<double> Energies;
    std::vector<double> Xsections;
} ;

// Cross-section-list is a collection of named Particle data
using XsecList = std::map<std::string,ParticleData>;

class A2ChannelManager
{
private:
    XsecList _XList;



public:
    using productlist_t = std::vector<ParticleTypeDatabase::Type>;


    A2ChannelManager(); // TODO: mask out decays {list}

    std::vector<ParticleTypeTreeDatabase::Channel> GetChannels() const;

    double TotalXsection(const double Egamma) const;

    double Xsection(const ParticleTypeTreeDatabase::Channel& channel, const double Egamma) const;

    std::string GetPlutoProductString(const ParticleTypeTreeDatabase::Channel& channel) const;

    struct beamTargetProducts_t{
        const std::string Beam;
        const std::string Target;
        const std::string Product;
        beamTargetProducts_t(const std::string& beam,
                             const std::string& target,
                             const std::string& products):
            Beam(beam), Target(target), Product(products){}
    };
    beamTargetProducts_t MakeBeamTargetProduct(const ParticleTypeTreeDatabase::Channel& channel ) const;


};

}
}
}
#endif 

