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
    using ChannelSelector_t = std::function<bool(const ParticleTypeTreeDatabase::Channel&)>;
    enum class Selection {
        All,
        gpBeamTarget,
        gnBeamTarget
    };
    static ChannelSelector_t GetSelector(const Selection& selection);


    static std::vector<ParticleTypeTreeDatabase::Channel> GetProductionChannels(const ChannelSelector_t& selector
                                                                                = GetSelector(Selection::All));
    static double TotalXsection(const double Egamma,
                                const ChannelSelector_t& selector = GetSelector(Selection::All));

    static double Xsection(const ParticleTypeTreeDatabase::Channel& channel,
                           const double Egamma);

    static std::string GetPlutoProductString(const ParticleTypeTreeDatabase::Channel& channel);
    static std::string GetPlutoBeamString(const ParticleTypeTreeDatabase::Channel& channel);
    static std::string GetPlutoTargetString(const ParticleTypeTreeDatabase::Channel& channel);



};

}
}
}
