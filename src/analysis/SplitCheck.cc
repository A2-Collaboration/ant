#include "SplitCheck.h"
#include "plot/Histogram.h"
#include <algorithm>
#include "TH1D.h"
#include "TH2D.h"
#include "plot/root_draw.h"
#include "utils/matcher.h"
#include "utils/combinatorics.h"

using namespace std;
using namespace ant;

ant::analysis::SplitCheck::SplitCheck()
{
    HistogramFactory::SetName("SplitCheck");
    HistogramFactory::SetLoopColors(true);

    const BinSettings gamma_energy(1000, 0.0, 1000.0);
    const BinSettings angle_bins(180, 0.0, 90.0);
    const BinSettings theta_bins(180, 0.0, 180.0);
    const BinSettings phi_bins(360, -180.0, 180.0);

    const int max=5;
    for(int gammas_event=2; gammas_event<=max; ++gammas_event) {
        gamma_rank[gammas_event].reserve(gammas_event);
        HistogramFactory::ResetColors();
        for(int gamma=0; gamma<gammas_event; ++gamma) {
            gamma_rank[gammas_event].push_back(
                        HistogramFactory::Make1D(to_string(gammas_event)+" gamma event: energy of gamma "+to_string(gamma),
                                  "E [MeV]",
                                  "",
                                  gamma_energy)
                        );
        }
    }

    big_small_angle = HistogramFactory::Make1D("Big Cluster/Small Cluster angle",
                                "angle [#circ]",
                                "",
                                angle_bins,
                                "big_small_angle");

    small_theta_phi = HistogramFactory::Make2D("Position of small clusters",
                                "#theta [#circ]",
                                "#phi [#circ]",theta_bins,phi_bins,"small_theta_phi");
    big_theta_phi = HistogramFactory::Make2D("Position of big clusters",
                                "#theta [#circ]",
                                "#phi [#circ]",theta_bins,phi_bins,"big_theta_phi");

    IMsmall_theta_phi = HistogramFactory::Make2D("Position of small clusters (IM selected)",
                                "#theta [#circ]",
                                "#phi [#circ]",theta_bins,phi_bins,"IMsmall_theta_phi");
    IMbig_theta_phi = HistogramFactory::Make2D("Position of big clusters (IM selected)",
                                "#theta [#circ]",
                                "#phi [#circ]",theta_bins,phi_bins,"IMbig_theta_phi");
}

void ant::analysis::SplitCheck::ProcessEvent(const ant::Event &event)
{
    if(event.Reconstructed().TriggerInfos().CBEenergySum()<550.0) return;

    ParticleList gammas = event.Reconstructed().Particles().Get(ParticleTypeDatabase::Photon);

    const size_t ngammas = gammas.size();


    if(ngammas==4) {
        ParticleList small_gammas;
        ParticleList large_gammas;

        for( auto& g :gammas ) {
            if(g->E() < 50.0 ) {
                small_gammas.push_back(g);
                small_theta_phi->Fill(g->Theta()*TMath::RadToDeg(), g->Phi()*TMath::RadToDeg());
            }
            else {
                large_gammas.push_back(g);
                big_theta_phi->Fill(g->Theta()*TMath::RadToDeg(), g->Phi()*TMath::RadToDeg());
            }
        }

        auto matches = utils::match1to1(large_gammas, small_gammas, [] ( const ParticlePtr& p1, const ParticlePtr& p2 ) {
                                            return p1->Angle(p2->Vect());
                                        });

        for(auto& big_small : matches) {
            big_small_angle->Fill(big_small.score*TMath::RadToDeg());
        }

    }

    for( auto c = makeCombination(gammas, 2); !c.Done(); ++c) {
        TLorentzVector m = *c.at(0) + *c.at(1);
        const double M = m.M();
        if(M<90.0) {
            IMsmall_theta_phi->Fill(c.at(0)->Theta()*TMath::RadToDeg(), c.at(0)->Phi()*TMath::RadToDeg());
            IMsmall_theta_phi->Fill(c.at(1)->Theta()*TMath::RadToDeg(), c.at(1)->Phi()*TMath::RadToDeg());
        } else {
            IMbig_theta_phi->Fill(c.at(0)->Theta()*TMath::RadToDeg(), c.at(0)->Phi()*TMath::RadToDeg());
            IMbig_theta_phi->Fill(c.at(1)->Theta()*TMath::RadToDeg(), c.at(1)->Phi()*TMath::RadToDeg());
        }
    }

    auto entry = gamma_rank.find(ngammas);
    if( entry == gamma_rank.end())
        return;

    auto& hists = entry->second;

    sort(gammas.begin(), gammas.end(), [] (const ParticlePtr& a, const ParticlePtr& b) { return a->E() > b->E();});

    for(size_t g=0;g<ngammas;++g) {
        hists.at(g)->Fill(gammas.at(g)->E());
    }





    /*
    // loop over the smaller clusters
    for(size_t g=ngammas_expected; g < gammas.size();++g) {
        const refRecParticle small = gammas.at(g);

        // loop over the bigger clusters
        for( size_t big_i=0; big_i<ngammas_expected; ++big_i) {

            //fill angle
            const double a = small->Angle(gammas.at(big_i)->Vect()) * TMath::RadToDeg();
            big_small_angle->Fill(a);
        }
    }
    */

}

void ant::analysis::SplitCheck::Finish()
{

}

void ant::analysis::SplitCheck::ShowResult()
{
    canvas c("gamma energy ranking");
    c << drawoption("nostack");

    for( auto& e : gamma_rank ) {
        const int ngammas = e.first;

        hstack s(to_string(ngammas) + "_gammas", to_string(ngammas)+"#gamma event: energy ranking");

        for(auto& h : e.second ) {
            s << h;
        }
        c << s;
    }
    c << drawoption("") << big_small_angle << endc;

    canvas c2 ("SplitCheck 2");
    c2 << drawoption("colz") << small_theta_phi << big_theta_phi << IMsmall_theta_phi << IMbig_theta_phi << endc;
}
