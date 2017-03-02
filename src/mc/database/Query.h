#pragma once

#include <vector>

#include "ProductionDataBase.h"


namespace ant
{
namespace mc
{
namespace data
{

// class with helper-functions for easier database acces
struct Query
{
    using ChannelSelector_t = std::function<bool(const ParticleTypeTreeDatabase::Channel&)>;

    /**
     * @brief Small database for often used selections of channels
     */
    enum class Selection {
        All,
        gpBeamTarget,
        gnBeamTarget
    };
    static ChannelSelector_t GetSelector(const Selection& selection);


    /**
     * @brief GetProductionChannels returns a vector with all channels which have cross-section data
     * @param selector : provide functions to contrain channels: channel -> bool
     * @return vector of TypeTree-channels
     */
    static std::vector<ParticleTypeTreeDatabase::Channel> GetProductionChannels(const ChannelSelector_t& selector
                                                                                = GetSelector(Selection::All));
    /**
     * @brief TotalXsection calculates the total cross-section (summing over database)
     * @param Egamma incident energy in MeV
     * @param selector function: channel -> bool
     * @return total cross section in mub
     */
    static double TotalXsection(const double Egamma,
                                const ChannelSelector_t& selector = GetSelector(Selection::All));
    /**
     * @brief Xsection get cross section for channel
     * @param channel ParticleTypeTree - channel
     * @param Egamma incident energy in MeV
     * @return xross section in mub
     */
    static double Xsection(const ParticleTypeTreeDatabase::Channel& channel,
                           const double Egamma);

    static std::string GetPlutoProductString(const ParticleTypeTreeDatabase::Channel& channel);
    /**
     * @brief GetPlutoBeamString xtracts the beam particle from given channel and checks if beam is known
     * @param channel
     * @return Pluto name for beam particle
     */
    static std::string GetPlutoBeamString(const ParticleTypeTreeDatabase::Channel& channel);
    /**
     * @brief GetPlutoTargetString xtracts the target particle from given channel and checks if it's known
     * @param channel
     * @return Pluto name for target particle
     */
    static std::string GetPlutoTargetString(const ParticleTypeTreeDatabase::Channel& channel);

};

}
}
}
