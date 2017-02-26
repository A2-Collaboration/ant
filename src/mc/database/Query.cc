#include "Query.h"

#include "base/Logger.h"
#include "analysis/utils/particle_tools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::mc::data;



vector<ParticleTypeTreeDatabase::Channel> Query::GetProductionChannels()
{
    vector<ParticleTypeTreeDatabase::Channel> channels;
    for (const auto& ch: ProductionDataBase::XSections) channels.emplace_back(ch.first);
    return channels;
}

double Query::TotalXsection(const double Egamma)
{
    double sigmatot = 0;

    for ( auto ch: ProductionDataBase::XSections)
        sigmatot += Xsection(ch.first, Egamma);

    return sigmatot;
}

double Query::Xsection(const ParticleTypeTreeDatabase::Channel &channel, const double Egamma)
{
    if (ProductionDataBase::XSections.count(channel))
        return ProductionDataBase::XSections.at(channel)(Egamma);

    LOG(WARNING) << "Production Channel " << utils::ParticleTools::GetDecayString(
                        ParticleTypeTreeDatabase::Get(channel)) << " not in Database!";

    return std::numeric_limits<double>::quiet_NaN();
}

string Query::GetPlutoProductString(const ParticleTypeTreeDatabase::Channel &channel)
{
    const auto tree = ParticleTypeTreeDatabase::Get(channel);
    return utils::ParticleTools::GetPlutoProduction(tree);
}

pair<string,string> getBeamTarget(const ParticleTypeTreeDatabase::Channel &channel)
{
    pair<string,string> beamTarget;
    const auto& root = ParticleTypeTreeDatabase::Get(channel)->Get();
    if ( root == ParticleTypeDatabase::BeamProton )
        beamTarget.first = ParticleTypeDatabase::Photon.PlutoName();
    else if ( root == ParticleTypeDatabase::BeamNeutron )
        beamTarget.first = ParticleTypeDatabase::Photon.PlutoName();
    else throw std::runtime_error("Beam particle unknown to database");

    if ( root == ParticleTypeDatabase::BeamProton )
        beamTarget.second = ParticleTypeDatabase::Proton.PlutoName();
    else if ( root == ParticleTypeDatabase::BeamNeutron )
        beamTarget.second = ParticleTypeDatabase::Neutron.PlutoName();
    else throw std::runtime_error("Target unknwon to database");

    return beamTarget;
}

string Query::GetPlutoBeamString(const ParticleTypeTreeDatabase::Channel &channel)
{
   return getBeamTarget(channel).first;
}

string Query::GetPlutoTargetString(const ParticleTypeTreeDatabase::Channel &channel)
{
    return getBeamTarget(channel).second;
}
