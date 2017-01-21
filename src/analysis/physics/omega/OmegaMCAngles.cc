#include "OmegaMCAngles.h"
#include "base/ParticleTypeTree.h"
#include "utils/particle_tools.h"
#include "analysis/plot/root_draw.h"


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::std_ext;





OmegaMCAngles::OmegaMCAngles(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    channelTree(ParticleTypeTreeDatabase::Get(ParticleTypeTreeDatabase::Channel::Pi0Eta_4g))
{
    const BinSettings bins(360,0,180);

    angle_pi0_gg  = HistFac.makeTH1D("Angle between #gamma from #pi^{0}", "angle [#circ]", "", bins, "angle_pi0");
    angle_eta_gg  = HistFac.makeTH1D("Angle between #gamma from eta",     "angle [#circ]", "", bins, "angle_eta");
    angle_all     = HistFac.makeTH1D("Angle between all #gamma",          "angle [#circ]", "", bins, "angle_all");
    angle_min_all = HistFac.makeTH1D("Minimal angle between all #gamma",  "angle [#circ]", "", bins, "angle_min_all");
}

OmegaMCAngles::~OmegaMCAngles()
{

}


template <typename T>
class MinTracker {
protected:
    T value = std::numeric_limits<double>::infinity();
public:

    MinTracker() {}

    void Set(const T& v) noexcept {
        if( v < value)
            value = v;
    }

    T Get() const noexcept {
        return value;
    }

    operator T() const noexcept {
        return value;
    }
};


void OmegaMCAngles::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& particletree = event.MCTrue().ParticleTree;

    if(!particletree || !particletree->IsEqual(channelTree, utils::ParticleTools::MatchByParticleName))
        return;


    TParticleList true_particles(4);

    particletree->Map_nodes([&true_particles] (const TParticleTree_t& p) {

            if(p->Get()->Type() == ParticleTypeDatabase::Pi0 && p->Daughters().size() == 2) {
                true_particles[0] = p->Daughters().front()->Get();
                true_particles[1] = p->Daughters().back()->Get();
            }

            if(p->Get()->Type() == ParticleTypeDatabase::Eta && p->Daughters().size() == 2) {
                true_particles[2] = p->Daughters().front()->Get();
                true_particles[3] = p->Daughters().back()->Get();
            }

    });

    angle_pi0_gg->Fill(radian_to_degree(true_particles.at(0)->Angle(*true_particles.at(1))));
    angle_eta_gg->Fill(radian_to_degree(true_particles.at(2)->Angle(*true_particles.at(3))));

    MinTracker<double> min_angle;
    for(auto i=true_particles.cbegin(); i!=true_particles.end(); ++i) {
        for(auto j=next(i); j!=true_particles.end(); ++j) {

            const auto a = (*i)->Angle(**j);

            angle_all->Fill(radian_to_degree(a));

            min_angle.Set(a);
        }
    }

    angle_min_all->Fill(radian_to_degree(min_angle.Get()));

}

void OmegaMCAngles::Finish()
{

}

void OmegaMCAngles::ShowResult()
{
    canvas(GetName()) << angle_all << angle_min_all << angle_eta_gg << angle_pi0_gg << endc;
}

AUTO_REGISTER_PHYSICS(OmegaMCAngles)
