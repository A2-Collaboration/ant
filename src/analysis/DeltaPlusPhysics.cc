#include "DeltaPlusPhysics.h"

#include "plot/Histogram.h"
#include "Event.h"
#include "Particle.h"
#include "ParticleType.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <iostream>

using namespace ant;
using namespace ant::analysis;
using namespace std;

DeltaPlusPhysics::DeltaPlusPhysics(const string &name):
    Physics(name),
    prompt("DeltaPlus_prompt"),
    random("DeltaPlus_random"),
    diff("DeltaPlus_diff"),
    pi0_cut(110,150),
    prompt_window(-8,8),
    random_window(-16,16),
    target(0,0,0,ParticleTypeDatabase::Proton.Mass())
{
    cout << "DeltaPlusPhysics:\n";
    cout << "Prompt window: " << prompt_window << " ns\n";
    cout << "Random window: " << random_window << " ns\n";
    cout << "Pi0 cut: " << pi0_cut << " MeV\n";
}

void DeltaPlusPhysics::ProcessEvent(const Event &event)
{
    ParticleList photons;
    ParticleList protons;

    for( auto& particle : event.Reconstructed().Particles().GetAll() ) {

        if( particle->Type() ==  ParticleTypeDatabase::Photon )
            photons.emplace_back(particle);
        else if ( particle->Type() == ParticleTypeDatabase::Proton )
            protons.emplace_back(particle);

    }

    for( auto& track : event.Reconstructed().Tracks() ) {
        prompt["pid"]->Fill(track->ClusterEnergy(), track->VetoEnergy());
    }

    for( auto& taggerhit : event.Reconstructed().TaggerHits()) {
        bool isPrompt = false;

        if( prompt_window.Contains(taggerhit->Time()) ) {
            isPrompt = true;
        } else if( random_window.Contains(taggerhit->Time())) {
            isPrompt = false;
        } else
            continue;

        Histogm& h = isPrompt ? prompt : random;

        // some basic histograms
        h["nPart"]->Fill(event.Reconstructed().Particles().GetAll().size());
        h["tag_energy"]->Fill(taggerhit->PhotonEnergy());
        h["tag_time"]->Fill(taggerhit->Time());


        if(photons.size() == 2) {
            const Particle pi0 ( ParticleTypeDatabase::Pi0, *photons.at(0) + *photons.at(1));
            h["2gIM"]->Fill(pi0.M());

            if( pi0_cut.Contains( pi0.M()) ) {
                h["pi0angle_noboost"]->Fill(pi0.Theta());

                // construct beam photon 4-vector and, using this, the delta restframe
                const TLorentzVector beam(0, 0, taggerhit->PhotonEnergy(), taggerhit->PhotonEnergy());
                const TLorentzVector delta_beam(beam + target);
                const TVector3 boost = -(delta_beam.BoostVector());

                // boost pi0
                TLorentzVector pi0_ = pi0;
                pi0_.Boost(boost);

                // missing mass plot (should peak at proton)
                TLorentzVector mmp = delta_beam - pi0;
                h["mmp"]->Fill( mmp.M() );

                // plot boosted pi0 angle, our desired plot!
                h["pi0angle"]->Fill( cos(pi0_.Theta()) );
                h["pi0angle_noboost"]->Fill( cos(pi0.Theta()));
                h["pi0angle_tagged"]->Fill( cos(pi0_.Theta()), beam.E());

                // have a look at proton reconstruction quality
                // by creating a delta
                if(protons.size()==1) {
                    TLorentzVector delta = pi0 + *protons.at(0);

                    TLorentzVector delta_ = delta;
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


void DeltaPlusPhysics::Histogm::AddHistogram(const string &name, const string &title, const string &x_label, const string &y_label, const int x_bins_n, const double x_bins_low, const double x_bins_up)
{

    // setup one dimensional histogram TH1D
    h_title[name] = title;

    h[name] = HistogramFactory::Default().Make1D(
                title,
                x_label,
                y_label,
                BinSettings(x_bins_n, x_bins_low, x_bins_up));

}

void DeltaPlusPhysics::Histogm::AddHistogram(const string &name, const string &title, const string &x_label, const string &y_label, const int x_bins_n, const double x_bins_low, const double x_bins_up, const int y_bins_n, const double y_bins_low, const double y_bins_up)
{

    // setup two dimensional histogram TH2D
    h_title[name] = title;
    h[name] = HistogramFactory::Default().Make2D(
                title,
                x_label,
                y_label,
                BinSettings(x_bins_n, x_bins_low, x_bins_up),
                BinSettings(y_bins_n, y_bins_low, y_bins_up));
}

DeltaPlusPhysics::Histogm::Histogm(const string &prefix)
{

    AddHistogram("nPart", "number of particles",
                   "number of particles / event", "",
                   10, 0, 10); // 10 bins from 0 to 10

    AddHistogram("pid", "PID Bananas",
                   "CB Energy [MeV]", "dE [MeV]",
                   100,0,450, // 100 bins from 0 to 450 in x
                   100,0,20   // 100 bins from 0 to 20  in y
                   );

    AddHistogram("2gIM", "2#gamma invariant mass",
                   "M_{#gamma #gamma} [MeV]", "",
                   100,0,300);

    AddHistogram("tag_time", "Tagger time",
                   "t [ns]", "",
                   100,-50,50);

    AddHistogram("tag_energy", "Tagged Photon Energy",
                   "E_{#gamma} [MeV]", "",
                   100,100,450);

    AddHistogram("mmp", "Missing Mass Proton",
                   "MM_{p} [MeV]", "",
                   100,600,1100);

    AddHistogram("pi0angle", "#pi^{0} #Theta angle (boosted)",
                   "cos(#theta_{#pi^{0}})", "",
                   180,-1,1);

    AddHistogram("pi0angle_noboost", "#pi^{0} #Theta angle (not boosted)",
                   "cos(#theta_{#pi^{0}})", "",
                   180,-1,1);

    AddHistogram("pi0angle_tagged", "#pi^{0} #Theta angle (boosted) vs tagged E",
                   "E_{#gamma} [MeV]", "cos(#theta_{pi^{0}})",
                   180,-1,1,14,110,300);

    AddHistogram("delta_pz", "#Delta^{+} momentum magnitude (boosted)",
                   "p_{#Delta^{+}} [MeV]","",
                   100,0,300);

    AddHistogram("delta_IM", "#Delta^{+} Invariant mass",
                   "M_{#Delta^{+}} [MeV]","",
                   100,800,1500);
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
