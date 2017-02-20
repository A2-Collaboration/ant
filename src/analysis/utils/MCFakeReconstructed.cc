#include "MCFakeReconstructed.h"

#include "utils/particle_tools.h"

#include "expconfig/detectors/CB.h"
#include "expconfig/detectors/PID.h"
#include "expconfig/detectors/TAPS.h"
#include "expconfig/detectors/TAPSVeto.h"

#include "base/Logger.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/misc.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::utils;

template<class T>
T find_detector() {
    using element_t = typename T::element_type;
    return ExpConfig::Setup::GetDetector<element_t>();
}

MCFakeReconstructed::MCFakeReconstructed(bool fakeComplete4Pi) :
    FakeComplete4Pi(fakeComplete4Pi),
    cb(find_detector<decltype(cb)>()),
    pid(find_detector<decltype(pid)>()),
    taps(find_detector<decltype(taps)>()),
    tapsveto(find_detector<decltype(tapsveto)>())
{
    LOG(WARNING) << "MCFakeReconstructed in use" << (FakeComplete4Pi ? ", with 4pi complete" : "");
}

MCFakeReconstructed::~MCFakeReconstructed()
{

}



unsigned find_closest_ch(const vec3& pos, const Detector_t& detector) {
    // assume to have at least one channel...
    unsigned min_ch = 0;
    double min_angle = pos.Angle(detector.GetPosition(min_ch));
    for(unsigned ch=1;ch<detector.GetNChannels();ch++) {
        const double angle = pos.Angle(detector.GetPosition(ch));
        if(angle < min_angle) {
            min_ch = ch;
            min_angle = angle;
        }
    }
    return min_ch;
}


void do_calo_veto(const TParticle& p,
                  const ClusterDetector_t& calo,
                  const Detector_t& veto,
                  TParticleList& list)
{
    auto calo_ch = find_closest_ch(p.p, calo);

    TClusterList clusters;

    clusters.emplace_back(
                p.p.Unit(),
                p.Ek(),
                0, // timing at zero for now
                calo.Type,
                calo_ch
                );
    Detector_t::Any_t type = calo.Type;


    /// \todo calculate better vetoE for charged particles?
    /// just assume 1 MeV for now...
    const double vetoE = p.Type().Charged() ? 1 : 0;

    if(vetoE>0) {
        auto veto_ch = find_closest_ch(p.p, veto);
        clusters.emplace_back(
                    veto.GetPosition(veto_ch).Unit(),
                    vetoE,
                    0, // timing zero for now
                    veto.Type,
                    veto_ch
                    );
        type |= veto.Type;
    }

    auto candidate = make_shared<TCandidate>(
                         type,
                         p.Ek(),
                         p.p.Theta(),
                         p.p.Phi(),
                         0,
                         0, // cluster size 0 for now...
                         vetoE,
                         std_ext::NaN, // no tracker energy
                         vetoE>0 ? TClusterList{std::prev(clusters.end(), 2), std::prev(clusters.end(), 1)}
                                 : TClusterList{std::prev(clusters.end())}
                                   );

    list.emplace_back(make_shared<TParticle>(p.Type(), candidate));
}

ParticleTypeList MCFakeReconstructed::Get(const TEventData& mctrue)
{
    TParticleList list;

    auto mctrue_particles = utils::ParticleTypeList::Make(mctrue.ParticleTree);

    for(const TParticlePtr& p : mctrue_particles.GetAll()) {
        auto type = geo.DetectorFromAngles(*p);
        if(type & Detector_t::Type_t::CB)
            do_calo_veto(*p, *cb, *pid, list);
        else if(type & Detector_t::Type_t::TAPS)
            do_calo_veto(*p, *taps, *tapsveto, list);
        else if(FakeComplete4Pi)
            do_calo_veto(*p, *cb, *pid, list); // assume lost particle is in CB

    }

    return ParticleTypeList::Make(list);
}