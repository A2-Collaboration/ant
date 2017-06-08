#include "pi0eta_lost_g.h"
#include "base/Logger.h"
#include "utils/ParticleTools.h"
#include <algorithm>
#include "base/std_ext/math.h"
#include "base/TH_ext.h"
#include "TTree.h"

#include "utils/Matcher.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;

Pi0EtaLostG::Pi0EtaLostG(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    lostPhtons(HistFac,"Photons"),
    lostProtons(HistFac,"Protons"),
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

    auto distance_function = [] (const TParticlePtr& p, const TCandidate& c) -> double {
        return p->Angle(c);
    };

    const auto getID = [] (const TParticle& p) {
        if(p.Type() == ParticleTypeDatabase::Photon)
            return 1;
        if(p.Type() == ParticleTypeDatabase::Proton)
            return 14;
        return -1;
    };

    t.lostV().clear();
    t.lostid().clear();
    for(const auto& m : utils::match2(mc_cands, re_cands, distance_function)) {
        const auto& p = *mc_cands.at(m.a);

        auto e = counters.find(addressof(p.Type()));
        if(e == counters.end()) {
            auto n = counters.emplace(make_pair(addressof(p.Type()), lostCounter(HistFac, p.Type().Name())));
            e = n.first;
        }
        lostCounter& c = e->second;

        c.FillProduced(p);

        if(!m.matched) {

            c.FillLost(p);

            t.lostV().emplace_back(p);
            t.lostid().emplace_back(getID(p));

        }
    }

    for(const auto& p : mc_cands) {
        if(p->Type() == ParticleTypeDatabase::Photon)
            lostPhtons.FillProduced(*p);
    }

    tree->Fill();

}

void Pi0EtaLostG::ShowResult()
{

    canvas c(GetName());
            c << TTree_drawable(tree,"lostV.Theta()*TMath::RadToDeg():lostV.Phi()*TMath::RadToDeg()","lostid==1")
            << TTree_drawable(tree,"lostV.Theta()*TMath::RadToDeg():lostV.Phi()*TMath::RadToDeg()","lostid==14")
            << drawoption("colz");
            for(auto& e : counters) {
                c << e.second.normalized;
            }
            c << endc;

}

void Pi0EtaLostG::Finish()
{
    for(auto& e : counters) {
        e.second.Finish(HistFac);
    }
}



Pi0EtaLostG::Tree_t::Tree_t()
{

}



Pi0EtaLostG::lostCounter::lostCounter(HistogramFactory &hf, const string& name):
    prefix(name)
{
    lostPhotons     = hf.makeTH2D("lost "+prefix,     "#theta [#circ]", "#phi[#circ]", BinSettings(180,0,180),BinSettings(360,-180,180), prefix+"_lost");
    producedPhotons = hf.makeTH2D("produced "+prefix, "#theta [#circ]", "#phi[#circ]", BinSettings(180,0,180),BinSettings(360,-180,180), prefix+"_produced");
}

void Pi0EtaLostG::lostCounter::Finish(HistogramFactory &hf)
{
    HistogramFactory::DirStackPush p(hf);
    normalized = TH_ext::Apply(lostPhotons, producedPhotons, [] (const double& a, const double& b) { if(b!=0.0) return a/b; return 0.0;}, prefix+"_norm");
    normalized->SetTitle((string("normalized ")+prefix).c_str());
}

void Pi0EtaLostG::lostCounter::FillLost(const TParticle &p)
{
    lostPhotons->Fill(radian_to_degree(p.Theta()), radian_to_degree(p.Phi()));
}

void Pi0EtaLostG::lostCounter::FillProduced(const TParticle &p)
{
    producedPhotons->Fill(radian_to_degree(p.Theta()), radian_to_degree(p.Phi()));
}

AUTO_REGISTER_PHYSICS(Pi0EtaLostG)
