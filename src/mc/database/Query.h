#pragma once

#include <vector>

#include "ProductionDataBase.h"


namespace ant
{
namespace mc
{
namespace data
{

    //TODO Define filters like: photon beam && proton target
struct Query
{
    static std::vector<ParticleTypeTreeDatabase::Channel> GetProductionChannels();

    static double TotalXsection(const double Egamma);
    static double Xsection(const ParticleTypeTreeDatabase::Channel& channel, const double Egamma);

    static std::string GetPlutoProductString(const ParticleTypeTreeDatabase::Channel& channel);
    static std::string GetPlutoBeamString(const ParticleTypeTreeDatabase::Channel& channel);
    static std::string GetPlutoTargetString(const ParticleTypeTreeDatabase::Channel& channel);



};

}
}
}
