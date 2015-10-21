#include "particle_tools.h"
#include "combinatorics.h"

#include "base/Logger.h"

#include "TH1.h"

#include <sstream>
#include <cassert>
#include <functional>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::data;

template<typename T>
string _GetDecayString(const shared_ptr<Tree<T>>& particletree, function<string(const T&)> to_string)
{
    if(!particletree)
        return "empty_unknown";

    stringstream s;

    // the head is the beam particle
    s << to_string(particletree->Get()) << " #rightarrow ";

    // ignore level==0 since it's the already handled beamparticle
    size_t lastlevel = 1;
    particletree->Map_level([&s, &lastlevel, to_string] (const T& p, size_t level) {
        if(level>0) {
            while(lastlevel<level) {
                s << "[ "; lastlevel++;
            }
            while(lastlevel>level) {
                s << "] "; lastlevel--;
            }
            assert(level == lastlevel);
            s << to_string(p) << " ";
        }
    });

    while(lastlevel-- > 1)
        s << "] ";

    return s.str();
}



string utils::ParticleTools::GetDecayString(const ParticleTree_t& particletree)
{
    return _GetDecayString<ParticlePtr>(particletree, [] (const ParticlePtr& p) { return p->Type().PrintName(); });
}

string utils::ParticleTools::GetDecayString(const ParticleTypeTree& particletypetree)
{
    return _GetDecayString<const ParticleTypeDatabase::Type&>(particletypetree, [] (const ParticleTypeDatabase::Type& t) { return t.PrintName(); });
}

string utils::ParticleTools::SanitizeDecayString(string decaystring)
{
    for(const auto c : {'(',')','[',']','{','}','^',' ','#'}) {
        std::replace( decaystring.begin(), decaystring.end(), c, '_');
    }
    std::replace(decaystring.begin(), decaystring.end(), '\'', 'p');
    std::replace(decaystring.begin(), decaystring.end(), '+', 'p');
    std::replace(decaystring.begin(), decaystring.end(), '-', 'm');
    return string("x")+decaystring;
}

string utils::ParticleTools::GetProductionChannelString(const data::ParticleTree_t& particletree)
{
    if(!particletree)
        return "empty_unknown";

    const auto& p = particletree->Get();

    stringstream s;

    s << p->Type().PrintName() << " #rightarrow";

    for(const auto& daughter : particletree->Daughters()) {
        s << " " << daughter->Get()->Type().PrintName();
    }

    return s.str();
}

const ParticlePtr utils::ParticleTools::FindParticle(const ant::ParticleTypeDatabase::Type& type, const ParticleList& particles)
{
    for(const auto& p : particles) {
        if(p->Type() == type) {
            return p;
        }
    }

    return nullptr;
}

void utils::ParticleTools::FillIMCombinations(TH1* h, unsigned n, const ParticleList& particles)
{
    for( auto comb = utils::makeCombination(particles,n); !comb.Done(); ++comb) {
         TLorentzVector sum(0,0,0,0);
         for(const auto& p : comb) {
             sum += *p;
         }
         h->Fill(sum.M());
    }
}

bool utils::ParticleTools::SortParticleByName(const data::ParticlePtr& a, const data::ParticlePtr& b)
{
    return a->Type().Name() < b->Type().Name();
}

bool utils::ParticleTools::MatchByParticleName(const data::ParticlePtr& a, const ant::ParticleTypeDatabase::Type& b)
{
    return a->Type().Name() == b.Name();
}
