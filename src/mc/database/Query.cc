#include "Query.h"

#include "base/Logger.h"
#include "analysis/utils/ParticleTools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::mc::data;



Query::ChannelSelector_t Query::GetSelector(const Query::Selection& selection)
{
    switch (selection) {
    case Query::Selection::gpBeamTarget:
        return [](const ParticleTypeTreeDatabase::Channel& ch)
            {
                return (ParticleTypeTreeDatabase::Get(ch)->Get()
                        == ParticleTypeDatabase::BeamProton);
            };
    case Query::Selection::gnBeamTarget:
        return [](const ParticleTypeTreeDatabase::Channel& ch)
            {
                return (ParticleTypeTreeDatabase::Get(ch)->Get()
                        == ParticleTypeDatabase::BeamNeutron);
            };
    case Query::Selection::All:
        return [](const ParticleTypeTreeDatabase::Channel&){return true;};
    }
    return [](const ParticleTypeTreeDatabase::Channel&){return true;};
}

vector<ParticleTypeTreeDatabase::Channel> Query::GetProductionChannels(const ChannelSelector_t &selector)
{
    vector<ParticleTypeTreeDatabase::Channel> channels;
    // higher order functions in c++ is a pain
    for_each(ProductionDataBase::XSections.begin(),
             ProductionDataBase::XSections.end(),
             [&selector,&channels](const ProductionDataBase::XSections_t::value_type& entry)
             {
                if (selector(entry.first))
                    channels.emplace_back(entry.first);
             } );
    return channels;
}

double Query::TotalXsection(const double Egamma,
                            const ChannelSelector_t& selector)
{
    double sigmatot = 0;

    for_each(ProductionDataBase::XSections.begin(),
             ProductionDataBase::XSections.end(),
             [&selector,&sigmatot,Egamma](const ProductionDataBase::XSections_t::value_type& entry)
             {
                if (selector(entry.first))
                    sigmatot += Xsection(entry.first, Egamma);
             } );

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
