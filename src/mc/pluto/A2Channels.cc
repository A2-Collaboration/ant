#include "A2Channels.h"
#include "mc/pluto/chgen.h"

#include <list>
#include <iostream>
#include <sstream>
#include <fstream>

#include "TCanvas.h"
#include "TH1D.h"

#include "base/Logger.h"
#include "analysis/utils/particle_tools.h"

// Pluto
#include "PDecayChannel.h"

using namespace std;
using namespace ant;
using namespace ant::mc::pluto;
using namespace ant::analysis;


A2ChannelManager::A2ChannelManager()
{
    _XList = LoadStdDecays();
}

double A2ChannelManager::Xsection(const ParticleTypeTreeDatabase::Channel& channel, const double Egamma) const
{
    if (ChannelDataBase::XSections.count(channel))
        return ChannelDataBase::XSections.at(channel)(Egamma);

    LOG(WARNING) << "Production Channel " << utils::ParticleTools::GetDecayString(
                        ParticleTypeTreeDatabase::Get(channel)) << " not in Database!";

    return std::numeric_limits<double>::quiet_NaN();
}

string A2ChannelManager::GetPlutoProductString(const ParticleTypeTreeDatabase::Channel& channel) const
{
     const auto tree = ParticleTypeTreeDatabase::Get(channel);

     return utils::ParticleTools::GetDecayString(tree,false);
}

A2ChannelManager::beamTargetProducts_t A2ChannelManager::MakeBeamTargetProduct(const ParticleTypeTreeDatabase::Channel& channel) const
{

}

double A2ChannelManager::TotalXsection(const double Egamma) const
{
    double sigmatot = 0;

    for ( auto ch: ChannelDataBase::XSections)
        sigmatot += Xsection(ch.first, Egamma);

    return sigmatot;
}






std::vector<ParticleTypeTreeDatabase::Channel> A2ChannelManager::GetChannels() const
{
    vector<ParticleTypeTreeDatabase::Channel> channels;
    for (const auto& ch: ChannelDataBase::XSections) channels.emplace_back(ch.first);
    return channels;
}

