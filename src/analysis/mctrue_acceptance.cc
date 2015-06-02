#include "analysis/mctrue_acceptance.h"
#include "plot/root_draw.h"

using namespace std;
using namespace ant;

analysis::MCTrueAcceptance::MCTrueAcceptance():
    Physics("McTrueAcceptance"), events_seen(0)
{
    detect = HistFac.makeHist<std::string>("Geometric Acceptance (MC True)","Particle Type","% events | all particles found", BinSettings(0),"detection");
    for( auto& type : ParticleTypeDatabase::DetectableTypes() ) {
        detect.Fill("all " + type->PrintName() + " detected", 1.0);
        detect.Fill("all " + type->PrintName() + " in CB", 1.0);
        detect.Fill("all " + type->PrintName() + " in TAPS", 1.0);
        detect.Fill("all " + type->PrintName() + " in CB & TAPS", 1.0);
    }
    detect.GetRootHistogram()->Reset();
}

analysis::MCTrueAcceptance::det_hit_count_t analysis::MCTrueAcceptance::AllAccepted(const ParticleList& particles) {
    det_hit_count_t acc;

    for( auto& p : particles ) {
        const detector_t d = geo.DetectorFromAngles(*p);

        if( d == detector_t::NaI) {
            acc.cb++;
        } else if( d & (detector_t::BaF2 | detector_t::PbWO4) ) {
            acc.taps++;
        }
    }

    return acc;
}

void ant::analysis::MCTrueAcceptance::ProcessEvent(const ant::Event &event)
{

    for( auto& type : ParticleTypeDatabase::DetectableTypes() ) {

        const ParticleList& particles = event.MCTrue().Particles().Get(*type);

        if(particles.size() > 0) {

            det_hit_count_t hit = AllAccepted(particles);

            if( particles.size() == hit.taps + hit.cb)
                detect.Fill("all " + type->PrintName() + " detected", 1.0);

            if( particles.size() == hit.cb )
                detect.Fill("all " + type->PrintName() + " in CB", 1.0);
            else if( particles.size() == hit.taps)
                detect.Fill("all " + type->PrintName() + " in TAPS", 1.0);
            else if( particles.size() == hit.taps + hit.cb)
                detect.Fill("all " + type->PrintName() + " in CB & TAPS", 1.0);
        }

    }

    events_seen++;

}

void ant::analysis::MCTrueAcceptance::Finish()
{
    detect.Scale(100.0/events_seen);
}

void ant::analysis::MCTrueAcceptance::ShowResult()
{
    canvas("MCTrueAcceptance") << detect << endc;
}

