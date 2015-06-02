#include "omega.h"
#include "Particle.h"
#include "Event.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "plot/root_draw.h"
#include "plot/Histogram.h"
#include "utils/combinatorics.h"
#include "TaggerHit.h"
#include <string>
#include <iostream>
#include "plot/SmartHist.h"

using namespace std;

ant::SmartHist1<const TLorentzVector&> ant::analysis::Omega::makeInvMassPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, BinSettings bins, const std::string& name) {
    return HistFac.makeHist<const TLorentzVector&>(
                [] (const TLorentzVector& p) { return p.M();},
                title,
                xlabel, ylabel, bins, name);
}

ant::SmartHist1< std::pair<const TLorentzVector&, const TLorentzVector&> > ant::analysis::Omega::makeAngleDiffPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, BinSettings bins, const std::string& name) {
    return HistFac.makeHist<std::pair<const TLorentzVector&, const TLorentzVector&>>(
            [] (std::pair<const TLorentzVector&, const TLorentzVector&> particles) {
                 return particles.first.Angle(particles.second.Vect())* TMath::RadToDeg();
            },
            title,
            xlabel, ylabel, bins, name);
}

ant::analysis::Omega::Omega(const string &name, bool mctrue, const ant::mev_t energy_scale):
    Physics(name),
    eta_im_cut(   IntervalD::CenterWidth( ParticleTypeDatabase::Eta.Mass(), 50.0)),
    pi0_im_cut( IntervalD::CenterWidth(ParticleTypeDatabase::Pi0.Mass(),20.0)),
    omega_im_cut( IntervalD::CenterWidth( ParticleTypeDatabase::Omega.Mass(), 80.0)),
    tagger_energy_cut(1420, 1575),
    target(0.0, 0.0, 0.0, ParticleTypeDatabase::Proton.Mass()),
    run_on_true(mctrue)
{
    const BinSettings energy_bins(1000, 0.0, energy_scale);
    const BinSettings p_MM_bins(1000, 500.0, 1500.0);
    const BinSettings angle_diff_bins(200,0.0,20.0);

    eta_IM      = makeInvMassPlot("2 #gamma IM (after omega cut)",  "M_{3#gamma}", "", energy_bins, "eta_IM");
    omega_IM    = makeInvMassPlot("3 #gamma IM (->#omega)",         "M_{3#gamma}", "", energy_bins, "omega_IM");
    p_MM        = makeInvMassPlot("MM",                             "MM [MeV]",    "", p_MM_bins,   "omega_MM");

    omega_rec_multi = HistFac.makeHist<int>(
                ParticleTypeDatabase::Omega.PrintName() + " Reconstruction Multiplicity",
                "n",
                "",
                BinSettings(5));

    nr_ngamma = HistFac.makeHist<int>(
                "Not reconstructed: number of photons",
                "number of photons/event",
                "",
                BinSettings(16));

    nr_2gim = HistFac.makeHist<const TLorentzVector&>([] (const TLorentzVector& v) { return v.M(); },
                "Not reconstructed: 2#gamma IM",
                "M_{2#gamma} [MeV]",
                "",
                energy_bins);

    nr_3gim = HistFac.makeHist<const TLorentzVector&>([] (const TLorentzVector& v) { return v.M(); },
                "Not reconstructed: 3#gamma IM",
                "M_{3#gamma} [MeV]",
                "",
                energy_bins);

    step_levels = HistFac.makeHist<std::string>(
                "Check pass count",
                "Check",
                "# passed",
                BinSettings(10));

    omega_mc_rec_angle = makeAngleDiffPlot(
                ParticleTypeDatabase::Omega.PrintName()+" MC/Rec angle",
                "angle [#circ]",
                "# / " + to_string(angle_diff_bins.BinWidth())+" #circ",
                angle_diff_bins,
                ParticleTypeDatabase::Eta.Name()+"_mc_rec_angle"
                );
    n=0;
}

template <class InputIterator, class T>
T sum (InputIterator first, InputIterator last, T init) {
    while (first!=last) {
        init += **first;
        ++first;
    }
    return std::move(init);
}

template <class C, class T>
T sum (const C& data, T init) {
    return std::move(sum(data.begin(), data.end(), init));
}

void ant::analysis::Omega::ProcessEvent(const ant::Event &event)
{

    const Event::Data& data = (run_on_true) ? event.MCTrue() : event.Reconstructed();

    step_levels.Fill("0 Events Seen");

 //   if(data.TriggerInfos().CBEenergySum()<550.0)
 //       return;

    step_levels.Fill("1 ESum Cut passed");

    const ParticleList& photons = data.Particles().Get(ParticleTypeDatabase::Photon);

    if(photons.size()<3)
        return;

    step_levels.Fill("2 NPhotons 3+");

    unsigned int n_omega_found = 0;

    for( auto comb = makeCombination(photons,3); !comb.Done(); ++comb) {

        ParticleList ggg;
        ggg.assign(comb.begin(),comb.end());

        TLorentzVector omega = *comb.at(0)+*comb.at(1)+*comb.at(2);

        if( omega_im_cut.Contains(omega.M())) {
            step_levels.Fill("3 #omega IM cut passed");

            for( auto gcomb = makeCombination(ggg,2); !gcomb.Done(); ++gcomb) {

                TLorentzVector g1(*gcomb.at(0));
                TLorentzVector g2(*gcomb.at(1));
                TLorentzVector eta = g1 + g2;

                eta_IM.Fill(eta);

                if(eta_im_cut.Contains(eta.M())) {

                    step_levels.Fill("4 #eta IM cut passed");
                    omega_IM.Fill(omega);
                    n_omega_found++;

                    for( auto& taggerhit : data.TaggerHits() ) {
                        if( tagger_energy_cut.Contains(taggerhit->PhotonEnergy())) {
                            TLorentzVector p = taggerhit->PhotonBeam() + target - omega;
                            p_MM.Fill(p);
                        }
                    }

                }
            }

        }

    }
    omega_rec_multi.Fill(n_omega_found);

    if(n_omega_found == 0) {
        nr_ngamma.Fill(photons.size());

        for( auto comb = makeCombination(photons,3); !comb.Done(); ++comb) {
            TLorentzVector m = sum(comb, TLorentzVector() );
            nr_3gim.Fill(m);
        }

        for( auto comb = makeCombination(photons,2); !comb.Done(); ++comb) {
            TLorentzVector m = sum(comb, TLorentzVector() );
            nr_2gim.Fill(m);
        }
    }

}


void ant::analysis::Omega::Finish()
{

}


void ant::analysis::Omega::ShowResult()
{
    canvas("Omega (Reconstructed)" + string((run_on_true ? " (True)" : "(Rec)"))) << omega_IM << eta_IM << p_MM << step_levels << omega_rec_multi << omega_mc_rec_angle << endc;
    canvas("Omega (Not Reconstructed)"+ string((run_on_true ? " (True)" : "(Rec)"))) << nr_ngamma << nr_2gim << nr_3gim << endc;

}



ant::SmartHist1<const TLorentzVector&> ant::analysis::Omega2::makeInvMassPlot(const std::string& title, const std::string& xlabel, const std::string& ylabel, BinSettings bins, const std::string& name) {
    return HistFac.makeHist<const TLorentzVector&>(
                [] (const TLorentzVector& p) { return p.M();},
                title,
                xlabel, ylabel, bins, name);
}


double ant::analysis::Omega2::calcEnergySum(const ParticleList &particles) const
{
    double esum = 0.0;

    for( const ant::ParticlePtr& track : particles) {
        if( geo.DetectorFromAngles(track->Theta(), track->Phi()) == detector_t::NaI ) {
            esum += track->Ek();
        }
    }

    return esum;
}

ant::ParticleList ant::analysis::Omega2::getGeoAccepted(const ant::ParticleList &p) const
{
    ParticleList list;
    for( auto& particle : p) {
        if( geo.DetectorFromAngles(particle->Theta(), particle->Phi()) != detector_t::None )
                list.emplace_back(particle);
    }
    return list;
}

ant::analysis::Omega2::Omega2(const string &name, bool mctrue):
    Physics(name),
    run_on_true(mctrue)
{
    settings.esum_threshold = 550.0;
    settings.omega_IM_cut = IntervalD(700, 850);
    settings.eta_IM_cut =   IntervalD(480, 620);
    settings.pi0_IM_cut =   IntervalD::CenterWidth( ParticleTypeDatabase::Pi0.Mass(), 40.0);

    step_levels = HistFac.makeHist<std::string>("Steps","","Pass count",BinSettings(0),"steps");
    const BinSettings energy_bins(1000);

    omega_IM    = makeInvMassPlot("#omega IM (after cut)",   "M_{#omega} [MeV]", "", energy_bins, "omega_IM");
    eta_IM      = makeInvMassPlot("#eta IM (after cut)",     "M_{#eta} [MeV]",   "", energy_bins, "eta_IM");
    p_MM        = makeInvMassPlot("MM",                      "MM [MeV]",         "", energy_bins, "MM");
    pi0_IM      = makeInvMassPlot("#pi^{0} IM (after cut)",  "M_{#pi^{0}} [MeV]", "", energy_bins, "pi0_IM");

    gggIM       = makeInvMassPlot("3 #gamma IM (before cut)",  "M_{3 #gamma} [MeV]", "", energy_bins, "ggg_IM");
     ggIM       = makeInvMassPlot("2 #gamma IM (before cut)",  "M_{2 #gamma} [MeV]", "", energy_bins,  "gg_IM");

    found_candidates = HistFac.makeHist<std::string>("Candidates found per Event","","",BinSettings(0),"candidates");

}

class BestParticleKeeper {
protected:
    ant::ParticlePtr particle;
    double current_dist;

    double dist(const ant::ParticlePtr& p) {
        return fabs(p->M() - p->Type().Mass());
    }

public:

    BestParticleKeeper(): current_dist(std::numeric_limits<double>::infinity()) {}

    void LookAt(ant::ParticlePtr& p) {
        double d = dist(p);
        if( !hasBest()) {
            particle = p;
            current_dist = d;
        } else {
            if(d < current_dist) {
                particle = p;
                current_dist = d;
            }
        }
    }

    ant::ParticlePtr GetBest() { return particle; }
    double GetBestDist() const { return current_dist; }
    bool hasBest() const { return (particle.get() != nullptr); }

};

void ant::analysis::Omega2::ProcessEvent(const ant::Event &event)
{
    bool omega_im_cut_passed(false);
    bool submeson_found(false);

    const Event::Data& data = (run_on_true) ? event.MCTrue() : event.Reconstructed();

    step_levels.Fill("Events Seen");

    if( data.Particles().Get(ParticleTypeDatabase::Photon).size() != 3 )
        return;

    step_levels.Fill("n_#gamma == 3");

    const double CBESum = run_on_true ? calcEnergySum(data.Particles().Get(ParticleTypeDatabase::Photon)) : data.TriggerInfos().CBEenergySum();

    if( CBESum < settings.esum_threshold )
        return;

    step_levels.Fill("CBESum");

    const ParticleList photons = getGeoAccepted(data.Particles().Get(ParticleTypeDatabase::Photon));
    const ParticleList protons = getGeoAccepted(data.Particles().Get(ParticleTypeDatabase::Proton));

    if(photons.size() != 3)
        return;

    step_levels.Fill("Geo Cut");

    vector<omega_decay> candidates;

    BestParticleKeeper best_eta;
    BestParticleKeeper best_pi0;

    for( auto comb = makeCombination(photons,3); !comb.Done(); ++comb) {

        ParticleList ggg;
        ggg.assign(comb.begin(),comb.end());

        TLorentzVector gggState = *comb.at(0)+*comb.at(1)+*comb.at(2);
        gggIM.Fill(gggState);

        if( settings.omega_IM_cut.Contains(gggState.M()) ) {
            omega_IM.Fill(gggState);

            if(!omega_im_cut_passed) { step_levels.Fill("3 #omega found"); omega_im_cut_passed=true; }

            for( auto gcomb = makeCombination(ggg,2); !gcomb.Done(); ++gcomb) {

                const TLorentzVector g1(*gcomb.at(0));
                const TLorentzVector g2(*gcomb.at(1));
                const TLorentzVector ggState = g1 + g2;

                ggIM.Fill(ggState);

                if(settings.eta_IM_cut.Contains(ggState.M())) {

                    candidates.emplace_back( omega_decay(
                                                 ParticlePtr(new Particle(ParticleTypeDatabase::Omega, gggState)),
                                                 ParticlePtr(new Particle(ParticleTypeDatabase::Eta,    ggState))));
                    if(!submeson_found) { step_levels.Fill("#eta/#pi^{0} found"); submeson_found=true; }

                    best_eta.LookAt( candidates.back().meson2 );

                } else if( settings.pi0_IM_cut.Contains(ggState.M())) {

                    candidates.emplace_back( omega_decay(
                                                 ParticlePtr(new Particle(ParticleTypeDatabase::Omega, gggState)),
                                                 ParticlePtr(new Particle(ParticleTypeDatabase::Pi0,    ggState))));
                    if(!submeson_found) { step_levels.Fill("#eta/#pi^{0} found"); submeson_found=true; }

                    best_pi0.LookAt( candidates.back().meson2 );

                }
            }
        }
    }

    //either eta or pi0, not both found
    if(best_eta.hasBest() ^ best_pi0.hasBest()) {
        if(best_eta.hasBest()) {
            eta_IM.Fill(*best_eta.GetBest().get());
        } else {
            pi0_IM.Fill(*best_pi0.GetBest().get());
        }
    }


   sort( candidates.begin(),
         candidates.end(),
         [] (const omega_decay& d1, const omega_decay& d2) { return d1.meson2->Type().Mass() > d2.meson2->Type().Mass(); });


    if( !candidates.empty() ) {
        string found_decays;
        for( const omega_decay& d : candidates) {
            found_decays += d.meson2->Type().PrintName();
        }
        found_candidates.Fill(found_decays);
    }

}


void ant::analysis::Omega2::Finish()
{

}

void ant::analysis::Omega2::ShowResult()
{
    canvas("Omega2")
            << gggIM << ggIM
            << omega_IM << eta_IM << pi0_IM
            << step_levels
            << found_candidates
            << endc;
}
