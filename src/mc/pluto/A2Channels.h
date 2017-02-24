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
    A2ChannelManager();

    std::vector<std::string> GetChannels() const;
    std::vector<ParticleTypeTreeDatabase::Channel> GetChannelsN() const;

    double Xsection(const std::string& name, const double Egamma) const;
    double TotalXsection(const double &Egamma) const;            // mask out decays {list}

    double Xsection(const ParticleTypeTreeDatabase::Channel& channel, const double Egamma) const;
    double TotalXsectionN(const double Egamma) const;
};

}
}
}
#endif 

