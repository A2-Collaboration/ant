#include "physics/common/ParticleIDCheck.h"

#include "base/std_ext/math.h"

#include <cmath>


using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace ant::analysis::data;

ParticleIDCheck::ParticleIDCheck(const std::string& name,PhysOptPtr opts):
    Physics(name, opts),
    mctrue(HistFac,"MCTrue"),
    rec(HistFac,"Rec")
{
    constexpr unsigned nPhiBins = 5;
    constexpr unsigned nThetaBins = 5;

    auto make_interval = [] (const double start, const double width) {
        return interval<double>(start,start+width);
    };

    for(unsigned i=0;i<nPhiBins;i++) {
        auto phi_range = make_interval(i*360.0/nPhiBins-180, 360.0/nPhiBins);
        for(unsigned j=0;j<nThetaBins;j++) {
            auto theta_range = make_interval(j*180.0/nThetaBins, 180.0/nThetaBins);

            bananas.emplace_back(phi_range, theta_range,
                                 HistFac.makeTH2D(std_ext::formatter()
                                                  << "Phi="<< phi_range
                                                  << "Theta="<< theta_range,
                                                  "Calorimeter Energy / MeV",
                                                  "Veto Energy / MeV",
                                                  BinSettings(400, 0, 800),
                                                  BinSettings(200, 0, 10)
                                                  )
                                 );
        }
    }
}

void ParticleIDCheck::ProcessEvent(const Event &event)
{
    mctrue.Fill(event.MCTrue());
    rec.Fill(event.Reconstructed());

    for(const auto& candidate : event.Reconstructed().Candidates()) {
        for(const auto& banana : bananas) {
            const auto& phi_range = std::get<0>(banana);
            const auto& theta_range = std::get<1>(banana);
            if(phi_range.Contains(std_ext::radian_to_degree(candidate->Phi())) &&
               theta_range.Contains(std_ext::radian_to_degree(candidate->Theta()))) {
                const auto& h = std::get<2>(banana);
                h->Fill(candidate->ClusterEnergy(), candidate->VetoEnergy());
            }
        }
    }
}


void ParticleIDCheck::Finish()
{

}


void ParticleIDCheck::ShowResult()
{
    canvas("ParticleIDCheck")
            <<  mctrue.hist << rec.hist
            << endc;
    canvas c_bananas("Veto Bananas");
    for(const auto& banana : bananas) {
        c_bananas << drawoption("colz") << std::get<2>(banana);
    }
    c_bananas << endc;
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

AUTO_REGISTER_PHYSICS(ParticleIDCheck)
