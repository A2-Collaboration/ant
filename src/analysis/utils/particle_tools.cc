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
using namespace ant::analysis::utils;

ParticleVars::ParticleVars(const TLorentzVector& lv, const ParticleTypeDatabase::Type& type) noexcept
{
    IM    = lv.M();
    Theta = std_ext::radian_to_degree(lv.Theta());
    Phi   = std_ext::radian_to_degree(lv.Phi());
    E     = lv.E() - type.Mass();
}

ParticleVars::ParticleVars(const TParticle& p) noexcept
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



string ParticleTools::GetDecayString(const TParticleTree_t& particletree)
{
    return _GetDecayString<TParticlePtr>(particletree, [] (const TParticlePtr& p) { return p->Type().PrintName(); });
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

string ParticleTools::GetProductionChannelString(const TParticleTree_t& particletree)
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

const TParticlePtr ParticleTools::FindParticle(const ant::ParticleTypeDatabase::Type& type, const TParticleList& particles)
{
    for(const auto& p : particles) {
        if(p->Type() == type) {
            return p;
        }
    }

    return nullptr;
}

const TParticlePtr ParticleTools::FindParticle(const ParticleTypeDatabase::Type& type,
                                              const TParticleTree_t& particletree,
                                              size_t maxlevel)
{
    if(!particletree)
        return nullptr;
    TParticlePtr p_found = nullptr;
    particletree->Map_level([&type, &p_found, maxlevel] (const TParticlePtr& p, size_t level) {
        if(level <= maxlevel && !p_found && p->Type() == type) {
            p_found = p;
        }
    });
    return p_found;
}

const TParticleList ParticleTools::FindParticles(const ParticleTypeDatabase::Type& type, const TParticleTree_t& particletree, size_t maxlevel)
{
    TParticleList list;
    if(!particletree)
        return list;
    particletree->Map_level([&type, &list, maxlevel] (const TParticlePtr& p, size_t level) {
        if(level <= maxlevel && p->Type() == type) {
            list.push_back(p);
        }
    });
    return list;
}

void ParticleTools::FillIMCombinations(TH1* h, unsigned n, const TParticleList& particles)
{
    FillIMCombinations([h] (double x) {h->Fill(x);}, n, particles);
}

void ParticleTools::FillIMCombinations(std::function<void(double)> filler, unsigned n, const TParticleList& particles)
{
    for( auto comb = makeCombination(particles,n); !comb.Done(); ++comb) {
         TLorentzVector sum(0,0,0,0);
         for(const auto& p : comb) {
             sum += *p;
         }
         filler(sum.M());
    }
}

bool ParticleTools::SortParticleByName(const TParticlePtr& a, const TParticlePtr& b)
{
    return a->Type().Name() < b->Type().Name();
}

bool ParticleTools::MatchByParticleName(const TParticlePtr& a, const ant::ParticleTypeDatabase::Type& b)
{
    return a->Type().Name() == b.Name();
}
