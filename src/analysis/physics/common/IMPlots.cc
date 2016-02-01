#include "IMPlots.h"
#include "utils/combinatorics.h"
#include "base/Logger.h"

#include "TH1D.h"
#include "TTree.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

IMPlots::IMPlots(const std::string& name, OptionsPtr opts):
    Physics(name, opts),
    m(8,{prs})
{
    prs.AddPromptRange({-2.5,2.5});
    prs.AddRandomRange({-15,-5});
    prs.AddRandomRange({  5,15});
    HistFac.SetTitlePrefix("aa");

    LOG(INFO) << "Promt Random Ratio = " << prs.Ratio();
    for(size_t i=0; i<m.size();++i) {
        m.at(i).MakeHistograms(HistFac,"IM_"+to_string(i+2),to_string(i+2)+" #gamma IM",BinSettings(1200),"IM [MeV]","");
    }
}

template <typename iter>
TLorentzVector sumlv(iter start, iter end) {
    TLorentzVector s;
    while(start != end ) {
        s += **(start);
        ++start;
    }
    return s;
}


void IMPlots::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& photons = event.Reconstructed->Particles.Get(ParticleTypeDatabase::Photon);

    for(unsigned n = MinNGamma(); n<MaxNGamma(); ++n) {
        for( auto comb = utils::makeCombination(photons,n); !comb.Done(); ++comb) {
            const TLorentzVector sum = sumlv(comb.begin(), comb.end());
                for(const auto& h : event.Reconstructed->Tagger.Hits) {
                    prs.SetTaggerHit(h.Time);
                    m.at(n - MinNGamma()).Fill(sum.M());
                }
            }
       }
}

void IMPlots::ShowResult()
{
    canvas c(GetName());

    for(auto& h : m) {

        c << h.subtracted;
    }

    c<< endc;
}


Symmetric2Gamma::Symmetric2Gamma(const string& name, OptionsPtr opts):
    Physics(name, opts)
{
    const BinSettings im(1600);
    h_symIM = HistFac.makeTH1D("2 #gamma IM, symmectic #gamma energies ("+to_string(perc*100)+"% Window)", "IM [MeV]", "", im, "symmetric2gamma");

    tree    = HistFac.makeTTree("tree");

    tree->Branch("IM",     &b_IM);
    tree->Branch("E",      &b_E);
    tree->Branch("E1",     &b_E1);
    tree->Branch("E2",     &b_E2);
    tree->Branch("Theta1", &b_theta1);
    tree->Branch("Phi1",   &b_phi1);
    tree->Branch("Theta2", &b_theta2);
    tree->Branch("Phi2",   &b_phi2);

}

Symmetric2Gamma::~Symmetric2Gamma()
{}

void Symmetric2Gamma::ProcessEvent(const TEvent& event, manager_t&)
{
    const auto& photons = event.Reconstructed->Particles.Get(ParticleTypeDatabase::Photon);
    for( auto comb = utils::makeCombination(photons,2); !comb.Done(); ++comb) {
        const TParticlePtr& g1 = comb.at(0);
        const TParticlePtr& g2 = comb.at(1);

        const auto Eavg = (g1->Ek()+ g2->Ek()) / 2.0;
        if(fabs(g1->Ek() - Eavg) < perc * Eavg) {
            const TLorentzVector sum = sumlv(comb.begin(), comb.end());

            b_IM = sum.M();
            b_E  = Eavg;

            b_E1     = g1->Ek();
            b_theta1 = g1->Theta();
            b_phi1   = g1->Phi();

            b_E2     = g2->Ek();
            b_theta2 = g2->Theta();
            b_phi2   = g2->Phi();

            tree->Fill();

            h_symIM->Fill(sum.M());
        }
    }

}

void Symmetric2Gamma::ShowResult()
{
    canvas(GetName())
            << h_symIM
            << endc;
}

AUTO_REGISTER_PHYSICS(IMPlots)
AUTO_REGISTER_PHYSICS(Symmetric2Gamma)
