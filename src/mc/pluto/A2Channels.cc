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


double A2ChannelManager::Xsection(const string &name, const double Egamma) const
{
    auto ch = _XList.find(name);
    if ( ch == _XList.end())
    {
        cerr << "  Warning:  '" << name << "' does not exist, returning 0 for Xsection." << endl;
        return 0;
    }

    // don't crash for unknown energies, just use last known value
    if (Egamma < ch->second.Energies.front()){
        return  ch->second.Xsections.front();
    }
    if (Egamma > ch->second.Energies.back()){
        return  ch->second.Xsections.back();
    }

    // interpolate:
    double xsec = ROOT::Math::Interpolator(ch->second.Energies,ch->second.Xsections).Eval(Egamma);

    //quickfix for Interpolator smoothing into negative numbers
    if (xsec < 0 ) xsec = 0;

    return xsec;
}

double A2ChannelManager::Xsection(const ParticleTypeTreeDatabase::Channel& channel, const double Egamma) const
{
    if (ChannelDataBase::XSections.count(channel))
        return ChannelDataBase::XSections.at(channel)(Egamma);

    LOG(WARNING) << "Production Channel " << utils::ParticleTools::GetDecayString(
                        ParticleTypeTreeDatabase::Get(channel)) << " not in Database!";

    return std::numeric_limits<double>::quiet_NaN();
}

double A2ChannelManager::TotalXsection(const double& Egamma) const
{
    double sigmatot = 0;

    for ( auto ch: _XList)
        sigmatot += Xsection(ch.first, Egamma);

    return sigmatot;
}

double A2ChannelManager::TotalXsectionN(const double Egamma) const
{
    double sigmatot = 0;

    for ( auto ch: ChannelDataBase::XSections)
        sigmatot += Xsection(ch.first, Egamma);

    return sigmatot;
}





A2ChannelManager::A2ChannelManager()
{

    _XList = LoadStdDecays();

}

vector<string> A2ChannelManager::GetChannels() const
{
    vector<string> names;
    for (const auto& ch: _XList) names.push_back(ch.first);
    return names;
}

std::vector<ParticleTypeTreeDatabase::Channel> A2ChannelManager::GetChannelsN() const
{
    vector<ParticleTypeTreeDatabase::Channel> channels;
    for (const auto& ch: ChannelDataBase::XSections) channels.emplace_back(ch.first);
    return channels;
}

