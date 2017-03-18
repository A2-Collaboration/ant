#include "DeltaPlusPhysics.h"

#include "utils/ParticleTools.h"
#include "base/ParticleType.h"

#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <iostream>

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;
using namespace std;

DeltaPlusPhysics::DeltaPlusPhysics(const string& name, OptionsPtr opts):
    Physics(name, opts),
    prompt(HistogramFactory("prompt", HistFac)),
    random(HistogramFactory("random", HistFac)),
    diff(HistogramFactory("diff", HistFac)),
    pi0_cut(110,150),
    prompt_window(-8,8),
    random_window(-16,16),
    target({0,0,0},ParticleTypeDatabase::Proton.Mass())
{
    cout << "DeltaPlusPhysics:\n";
    cout << "Prompt window: " << prompt_window << " ns\n";
    cout << "Random window: " << random_window << " ns\n";
    cout << "Pi0 cut: " << pi0_cut << " MeV\n";
}

void DeltaPlusPhysics::ProcessEvent(const TEvent& event, manager_t&)
{
    TParticleList photons;
    TParticleList protons;

    auto recon_particles = utils::ParticleTypeList::Make(event.Reconstructed().Candidates);

    for( auto& particle : recon_particles.GetAll() ) {

        if( particle->Type() ==  ParticleTypeDatabase::Photon )
            photons.emplace_back(particle);
        else if ( particle->Type() == ParticleTypeDatabase::Proton )
            protons.emplace_back(particle);

    }

    for(const auto& cand : event.Reconstructed().Candidates ) {
        prompt["pid"]->Fill(cand.CaloEnergy, cand.VetoEnergy);
    }

    for(const auto& taggerhit : event.Reconstructed().TaggerHits) {
        bool isPrompt = false;

        if( prompt_window.Contains(taggerhit.Time) ) {
            isPrompt = true;
        } else if( random_window.Contains(taggerhit.Time)) {
            isPrompt = false;
        } else
            continue;

        Histogm& h = isPrompt ? prompt : random;

        // some basic histograms
        h["nPart"]->Fill(recon_particles.GetAll().size());
        h["tag_energy"]->Fill(taggerhit.PhotonEnergy);
        h["tag_time"]->Fill(taggerhit.Time);


        if(photons.size() == 2) {
            const TParticle pi0 ( ParticleTypeDatabase::Pi0, *photons.at(0) + *photons.at(1));
            h["2gIM"]->Fill(pi0.M());

            if( pi0_cut.Contains( pi0.M()) ) {
                h["pi0angle_noboost"]->Fill(pi0.Theta());

                // construct beam photon 4-vector and, using this, the delta restframe
                const LorentzVec beam({0, 0, taggerhit.PhotonEnergy}, taggerhit.PhotonEnergy);
                const LorentzVec delta_beam(beam + target);
                const auto boost = -(delta_beam.BoostVector());

                // boost pi0
                LorentzVec pi0_ = pi0;
                pi0_.Boost(boost);

                // missing mass plot (should peak at proton)
                LorentzVec mmp = delta_beam - pi0;
                h["mmp"]->Fill( mmp.M() );

                // plot boosted pi0 angle, our desired plot!
                h["pi0angle"]->Fill( cos(pi0_.Theta()) );
                h["pi0angle_noboost"]->Fill( cos(pi0.Theta()));
                h["pi0angle_tagged"]->Fill( cos(pi0_.Theta()), beam.E);

                // have a look at proton reconstruction quality
                // by creating a delta
                if(protons.size()==1) {
                    LorentzVec delta = pi0 + *protons.at(0);

                    LorentzVec delta_ = delta;
                    delta_.Boost(boost);

                    h["delta_pz"]->Fill(delta_.P());
                    h["delta_IM"]->Fill(delta_.M());
                }
            }
        }

    }
}

void DeltaPlusPhysics::Finish()
{
    diff = prompt;
    diff.AddScaled(random, -1.0);
}

void DeltaPlusPhysics::ShowResult()
{
    diff.Draw();

}

DeltaPlusPhysics::Histogm::Histogm(HistogramFactory HistFac)
{

    auto insert_hist = [this] (TH1* hist) {
        h[hist->GetName()] = hist;
    };

    insert_hist(HistFac.makeTH1D("number of particles",
                                 "number of particles / event", "",
                                 BinSettings(10, 0, 10), "nPart")); // 10 bins from 0 to 10

    insert_hist(HistFac.makeTH2D("PID Bananas",
                                 "CB Energy [MeV]", "dE [MeV]",
                                 BinSettings(100,0,450),     // 100 bins from 0 to 450 in x
                                 BinSettings(100,0,20), // 100 bins from 0 to 20  in y
                                 "pid"));

    insert_hist(HistFac.makeTH1D("2#gamma invariant mass",
                                 "M_{#gamma #gamma} [MeV]", "",
                                 BinSettings(100,0,300), "2gIM"));

    insert_hist(HistFac.makeTH1D("Tagger time",
                                 "t [ns]", "",
                                 BinSettings(100,-50,50),"tag_time"));

    insert_hist(HistFac.makeTH1D("Tagged Photon Energy",
                                 "E_{#gamma} [MeV]", "",
                                 BinSettings(100,100,450),"tag_energy"));

    insert_hist(HistFac.makeTH1D("Missing Mass Proton",
                                 "MM_{p} [MeV]", "",
                                 BinSettings(100,600,1100),"mmp"));

    insert_hist(HistFac.makeTH1D("#pi^{0} #Theta angle (boosted)",
                                 "cos(#theta_{#pi^{0}})", "",
                                 BinSettings(180,-1,1),"pi0angle"));

    insert_hist(HistFac.makeTH1D("#pi^{0} #Theta angle (not boosted)",
                                 "cos(#theta_{#pi^{0}})", "",
                                 BinSettings(180,-1,1),"pi0angle_noboost"));

    insert_hist(HistFac.makeTH2D("#pi^{0} #Theta angle (boosted) vs tagged E",
                                 "E_{#gamma} [MeV]", "cos(#theta_{pi^{0}})",
                                 BinSettings(180,-1,1),
                                 BinSettings(14,110,300),"pi0angle_tagged"));

    insert_hist(HistFac.makeTH1D("#Delta^{+} momentum magnitude (boosted)",
                                 "p_{#Delta^{+}} [MeV]","",
                                 BinSettings(100,0,300),"delta_pz"));

    insert_hist(HistFac.makeTH1D("#Delta^{+} Invariant mass",
                                 "M_{#Delta^{+}} [MeV]","",
                                 BinSettings(100,800,1500),"delta_IM"));
}

void DeltaPlusPhysics::Histogm::Draw()
{
    {
        TCanvas* c = new TCanvas(Form("%s_c",pref.c_str()),pref.c_str());
        const int cols = ceil(sqrt(h.size()));
        const int rows = ceil((double)h.size()/(double)cols);
        c->Divide(cols,rows);
        int pad=1;
        for(auto i=h.begin(); i!=h.end(); ++i) {
            c->cd(pad);
            TH2D* h2 = dynamic_cast<TH2D*>(i->second);
            if(h2 != nullptr) {
                h2->Draw("colz");
            }
            else {
                i->second->Draw();
            }
            pad++;
        }
    }
}

DeltaPlusPhysics::Histogm &DeltaPlusPhysics::Histogm::operator*=(const Double_t factor)
{

    for(auto i=h.begin(); i!=h.end(); ++i) {
        i->second->Scale(factor);
    }
    return *this;

}

DeltaPlusPhysics::Histogm DeltaPlusPhysics::Histogm::operator=(const DeltaPlusPhysics::Histogm &other)
{

    for(auto i=h.begin(); i!=h.end(); ++i) {
        TH1* h = i->second;
        h->Reset();
        h->Add(other[i->first]);
    }
    return *this;

}
void DeltaPlusPhysics::Histogm::AddScaled(const DeltaPlusPhysics::Histogm &h2, const Double_t f)
{
    for(auto i=h.begin(); i!=h.end(); ++i) {
        i->second->Add(h2.h.at(i->first),f);
    }
}

AUTO_REGISTER_PHYSICS(DeltaPlusPhysics)
