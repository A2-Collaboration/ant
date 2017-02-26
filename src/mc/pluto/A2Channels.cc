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


ProductionTools::ProductionTools()
{
    _XList = LoadStdDecays();
}

double ProductionTools::Xsection(const ParticleTypeTreeDatabase::Channel& channel, const double Egamma) const
{
    if (ProductionDataBase::XSections.count(channel))
        return ProductionDataBase::XSections.at(channel)(Egamma);

    LOG(WARNING) << "Production Channel " << utils::ParticleTools::GetDecayString(
                        ParticleTypeTreeDatabase::Get(channel)) << " not in Database!";

    return std::numeric_limits<double>::quiet_NaN();
}

string ProductionTools::GetPlutoProductString(const ParticleTypeTreeDatabase::Channel& channel)
{
     const auto tree = ParticleTypeTreeDatabase::Get(channel);
     return utils::ParticleTools::GetPlutoProduction(tree);
}

string ProductionTools::GetBeam(const ParticleTypeTreeDatabase::Channel& channel)
{
    if ( ParticleTypeTreeDatabase::Get(channel)->Get() == ParticleTypeDatabase::BeamProton )
        return ParticleTypeDatabase::Photon.PlutoName();
    if ( ParticleTypeTreeDatabase::Get(channel)->Get() == ParticleTypeDatabase::BeamNeutron )
        return ParticleTypeDatabase::Photon.PlutoName();


    throw std::runtime_error("Unknown Beam");
    return "";
}

string ProductionTools::GetTarget(const ParticleTypeTreeDatabase::Channel& channel)
{
    if ( ParticleTypeTreeDatabase::Get(channel)->Get() == ParticleTypeDatabase::BeamProton )
        return ParticleTypeDatabase::Proton.PlutoName();
    if ( ParticleTypeTreeDatabase::Get(channel)->Get() == ParticleTypeDatabase::BeamNeutron )
        return ParticleTypeDatabase::Neutron.PlutoName();


    throw std::runtime_error("Unknown Beam");
    return "";
}





double ProductionTools::TotalXsection(const double Egamma) const
{
    double sigmatot = 0;

    for ( auto ch: ProductionDataBase::XSections)
        sigmatot += Xsection(ch.first, Egamma);

    return sigmatot;
}






std::vector<ParticleTypeTreeDatabase::Channel> ProductionTools::GetChannels() const
{
    vector<ParticleTypeTreeDatabase::Channel> channels;
    for (const auto& ch: ProductionDataBase::XSections) channels.emplace_back(ch.first);
    return channels;
}

