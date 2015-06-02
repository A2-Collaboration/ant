#include "RecoCheck.h"

#include "TH1D.h"
#include "utils/matcher.h"
#include "TLorentzVector.h"
#include "Particle.h"
#include "plot/root_draw.h"
#include "plot/Histogram.h"


using namespace std;
using namespace ant;



analysis::RecoCheck::RecoCheck():
    cb_angle(20.0*TMath::DegToRad(), 160.0*TMath::DegToRad())
{
    HistogramFactory::SetName("RecoCheck");
    const BinSettings angle_bins(300,0.0,30.0);
    const BinSettings npart_bins(10);

    angle_diff = HistogramFactory::Make1D( "MC/Rec Angle IM",
                           "angle [#circ]",
                           "# / " + to_string(angle_bins.BinWidth())+" #circ",
                           angle_bins,
                           "angle_diff"
                            );
    n_unmatched = HistogramFactory::Make1D("Unmatched particles / event",
                            "# unmatched / event","",npart_bins,"n_unmatched");
}

void analysis::RecoCheck::ProcessEvent(const Event &event)
{
    const ParticleList& mc = event.MCTrue().Particles().GetAll();
    const ParticleList& rec = event.Reconstructed().Particles().GetAll();

    ParticleList mc_in_cb;

    for( auto& p : mc ) {
        if( cb_angle.Contains(p->Theta()))
            mc_in_cb.push_back(p);
    }

    // find mc-reco matches
    auto matched = utils::match1to1(mc_in_cb, rec, [] ( const ParticlePtr& p1, const ParticlePtr& p2 ) {
                                        return p1->Angle(p2->Vect());
                                    });

    // cacluate angles for all matches
    for( auto& match : matched ) {
        const double a = match.score * TMath::RadToDeg();
        angle_diff->Fill(a);
    }

    // fill number of unmatched mc particles
    n_unmatched->Fill(mc_in_cb.size() - matched.size());

}

void analysis::RecoCheck::Finish()
{

}

void analysis::RecoCheck::ShowResult()
{
    canvas("RecoCheck") << angle_diff << n_unmatched << endc >> "A.pdf";
}

