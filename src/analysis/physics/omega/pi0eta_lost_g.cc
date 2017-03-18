#include "pi0eta_lost_g.h"
#include "base/Logger.h"
#include "utils/particle_tools.h"
#include <algorithm>
#include "base/std_ext/math.h"

#include "TTree.h"

#include "utils/Matcher.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;

Pi0EtaLostG::Pi0EtaLostG(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    tree(HistFac.makeTTree("tree"))
{
    t.CreateBranches(tree);
}

Pi0EtaLostG::~Pi0EtaLostG()
{
}

struct matchpair {
    size_t a;
    size_t b;
    double dist;

    matchpair(size_t index_a, size_t index_b, double distance) noexcept:
        a(index_a),
        b(index_b),
        dist(distance) {}

    bool operator> (const matchpair& o) const noexcept {
        return dist > o.dist;
    }

    bool operator< (const matchpair& o) const noexcept {
        return dist < o.dist;
    }
};

void Pi0EtaLostG::ProcessEvent(const TEvent& event, manager_t&)
{
    auto mctrue_particles = utils::ParticleTypeList::Make(event.MCTrue().ParticleTree);

    const auto& mc_cands = mctrue_particles.GetAll();
    const auto& re_cands = event.Reconstructed().Candidates;

    if(mc_cands.size() != 5 || re_cands.size() != 4)
        return;

    auto distance_function = [] (const TParticlePtr& p, const TCandidate& c) -> double {
        return p->Angle(c);
    };

    const auto matched = utils::match2(mc_cands, re_cands, distance_function);

    auto findUnmatched = [] (const decltype (matched)& l) -> const utils::matchpair& {
        for(const auto& m : l) {
            if(!m.matched) {
                return m;
            }
        }
        throw std::runtime_error("no unmatched particle found!");
    };

    const auto iUnmatched = findUnmatched(matched);


    const auto Unmatched = mc_cands.at(iUnmatched.a);

    t.E        = Unmatched->Ek();
    t.theta    = Unmatched->Theta();
    t.phi      = Unmatched->Phi();
    t.MinAngle = iUnmatched.dist;

    tree->Fill();

}



Pi0EtaLostG::Tree_t::Tree_t()
{

}

AUTO_REGISTER_PHYSICS(Pi0EtaLostG)
