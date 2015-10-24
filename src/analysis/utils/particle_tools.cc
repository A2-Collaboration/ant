#include "particle_tools.h"
#include "combinatorics.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"

#include "TH1.h"
#include "TTree.h"

#include <sstream>
#include <cassert>
#include <functional>
#include <limits>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::data;
using namespace ant::analysis::utils;

ParticleVars::ParticleVars(const TLorentzVector& lv, const ParticleTypeDatabase::Type& type) noexcept
{
    IM    = lv.M();
    Theta = std_ext::radian_to_degree(lv.Theta());
    Phi   = std_ext::radian_to_degree(lv.Phi());
    E     = lv.E() - type.Mass();
}

ParticleVars::ParticleVars(const Particle& p) noexcept
{
    IM    = p.M();
    Theta = std_ext::radian_to_degree(p.Theta());
    Phi   = std_ext::radian_to_degree(p.Phi());
    E     = p.Ek();
}

void ParticleVars::SetBranches(TTree* tree, const string& prefix)
{
    tree->Branch((prefix+"_IM").c_str(), &IM);
    tree->Branch((prefix+"_Theta").c_str(), &Theta);
    tree->Branch((prefix+"_Phi").c_str(),&Phi);
    tree->Branch((prefix+"_E").c_str(),  &E);
}

void ParticleVars::Clear()
{
    IM = numeric_limits<double>::quiet_NaN();
    Theta = numeric_limits<double>::quiet_NaN();
    Phi = numeric_limits<double>::quiet_NaN();
    E = numeric_limits<double>::quiet_NaN();
}


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



string ParticleTools::GetDecayString(const ParticleTree_t& particletree)
{
    return _GetDecayString<ParticlePtr>(particletree, [] (const ParticlePtr& p) { return p->Type().PrintName(); });
}

string ParticleTools::GetDecayString(const ParticleTypeTree& particletypetree)
{
    return _GetDecayString<const ParticleTypeDatabase::Type&>(particletypetree, [] (const ParticleTypeDatabase::Type& t) { return t.PrintName(); });
}

string ParticleTools::SanitizeDecayString(string decaystring)
{
    for(const auto c : {'(',')','[',']','{','}','^',' ','#'}) {
        std::replace( decaystring.begin(), decaystring.end(), c, '_');
    }
    std::replace(decaystring.begin(), decaystring.end(), '\'', 'p');
    std::replace(decaystring.begin(), decaystring.end(), '+', 'p');
    std::replace(decaystring.begin(), decaystring.end(), '-', 'm');
    return string("x")+decaystring;
}

string ParticleTools::GetProductionChannelString(const data::ParticleTree_t& particletree)
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

const ParticlePtr ParticleTools::FindParticle(const ant::ParticleTypeDatabase::Type& type, const ParticleList& particles)
{
    for(const auto& p : particles) {
        if(p->Type() == type) {
            return p;
        }
    }

    return nullptr;
}

const ParticlePtr ParticleTools::FindParticle(const ParticleTypeDatabase::Type& type,
                                              const ParticleTree_t& particletree,
                                              size_t maxlevel)
{
    if(!particletree)
        return nullptr;
    ParticlePtr p_found = nullptr;
    particletree->Map_level([&type, &p_found, maxlevel] (const ParticlePtr& p, size_t level) {
        if(level <= maxlevel && !p_found && p->Type() == type) {
            p_found = p;
        }
    });
    return p_found;
}

void ParticleTools::FillIMCombinations(TH1* h, unsigned n, const ParticleList& particles)
{
    for( auto comb = makeCombination(particles,n); !comb.Done(); ++comb) {
         TLorentzVector sum(0,0,0,0);
         for(const auto& p : comb) {
             sum += *p;
         }
         h->Fill(sum.M());
    }
}

bool ParticleTools::SortParticleByName(const data::ParticlePtr& a, const data::ParticlePtr& b)
{
    return a->Type().Name() < b->Type().Name();
}

bool ParticleTools::MatchByParticleName(const data::ParticlePtr& a, const ant::ParticleTypeDatabase::Type& b)
{
    return a->Type().Name() == b.Name();
}
