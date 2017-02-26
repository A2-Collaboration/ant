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

class ProductionTools
{
private:
    XsecList _XList;



public:
    using productlist_t = std::vector<ParticleTypeDatabase::Type>;


    ProductionTools(); // TODO: mask out decays {list}

    std::vector<ParticleTypeTreeDatabase::Channel> GetChannels() const;

    double TotalXsection(const double Egamma) const;

    double Xsection(const ParticleTypeTreeDatabase::Channel& channel, const double Egamma) const;

    static std::string GetPlutoProductString(const ParticleTypeTreeDatabase::Channel& channel);
    static std::string GetBeam(const ParticleTypeTreeDatabase::Channel& channel);
    static std::string GetTarget(const ParticleTypeTreeDatabase::Channel& channel);


};

}
}
}
#endif 

