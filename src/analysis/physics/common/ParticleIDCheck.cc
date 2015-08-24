#include "physics/common/ParticleIDCheck.h"
#include "plot/root_draw.h"
#include "data/Particle.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include <cmath>

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

ParticleIDCheck::ParticleIDCheck(PhysOptPtr opts):
    Physics("ParticleIDCheck", opts),
    mctrue(HistFac,"MCTrue"),
    rec(HistFac,"Rec")
{

}

void ParticleIDCheck::ProcessEvent(const Event &event)
{
    mctrue.Fill(event.MCTrue());
    rec.Fill(event.Reconstructed());
}


void ParticleIDCheck::Finish()
{

}


void ParticleIDCheck::ShowResult()
{
    canvas("ParticleIDCheck")
            <<  mctrue.hist << rec.hist
            << endc;
}



ParticleIDCheck::branch_hists::branch_hists(SmartHistFactory& HistFac, const string& name)
{
    hist = HistFac.makeTH1D(name+": particles in CB","","",BinSettings(10),name+"_particles");
    for(auto& pt : ParticleTypeDatabase::DetectableTypes()) {
        hist->Fill(pt->PrintName().c_str(),0);
    }
}

void ParticleIDCheck::branch_hists::Fill(const Event::Data& data)
{
    hist->Fill("unID", max(0,int(data.Candidates().size()) - int(data.Particles().GetAll().size())));

    for(auto& p: data.Particles().GetAll()) {
        hist->Fill(p->Type().PrintName().c_str(),1);
    }
}

AUTO_REGISTER_PHYSICS(ParticleIDCheck, "ParticleIDCheck")
