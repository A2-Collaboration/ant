#include "particle_tools.h"
#include "combinatorics.h"

#include "base/Logger.h"

#include "TH1.h"

#include <sstream>

using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace  std;

string utils::ParticleTools::GetDecayString(const ParticleTree_t& particletree)
{
    stringstream s;

    // the head is the beam particle
    s << particletree->get()->Type().PrintName() << " #rightarrow ";

    // ignore level==0 since its the already handled beamparticle
    size_t lastlevel = 1;
    particletree->maplevel([&s, &lastlevel] (const shared_ptr<Particle> p, size_t level) {
        if(level>0) {
            if(lastlevel < level)
                s << "[ ";
            else if(lastlevel > level)
                s << "] ";
            s << p->Type().PrintName() << " ";
            lastlevel = level;
        }
    });

    while(lastlevel-- > 1)
        s << "] ";

    return s.str();
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
    const auto& p = particletree->get();

    stringstream s;

    s << p->Type().PrintName() << " #rightarrow";

    for(const auto& daughter : particletree->Daughters()) {
        s << " " << daughter->get()->Type().PrintName();
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
