#include "particle_tools.h"
#include "combinatorics.h"

#include "TH1.h"

#include <sstream>

using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace  std;

string utils::ParticleTools::GetDecayString(const ParticleList& particles)
{
    stringstream s;
    if(! particles.empty()) {
        const auto& start = particles.front();
        Particle::RecPrint(start, s);
    }

    return s.str();
}

string utils::ParticleTools::SanitizeDecayString(string decaystring)
{
    for(const auto c : {'(',')','[',']','{','}','^',' ','#'}) {
        std::replace( decaystring.begin(), decaystring.end(), c, '_');
    }
    std::replace(decaystring.begin(), decaystring.end(), '\'', 'p');
    return string("x")+decaystring;
}

string utils::ParticleTools::GetProductionChannelString(const ParticleList& particles)
{
    const auto p = FindParticle(ParticleTypeDatabase::BeamTarget, particles);

    stringstream s;

    if(p) {

        s << p->Type().PrintName() << " #rightarrow";

        for(const auto& daughter : p->Daughters()) {
            s << " " << daughter->Type().PrintName();
        }

    } else {
        s << "???";
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