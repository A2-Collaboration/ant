#include "MCFakeReconstructed.h"

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

double get_avg_opening_angle(const ClusterDetector_t& detector) {
    double sum_angle = 0;
    for(unsigned ch=0;ch<detector.GetNChannels();ch++) {
        const ClusterDetector_t::Element_t* elem = detector.GetClusterElement(ch);
        double sum_angle_neigh = 0;
        for(unsigned ch_neigh : elem->Neighbours) {
            const ClusterDetector_t::Element_t* neigh = detector.GetClusterElement(ch_neigh);
            sum_angle_neigh += elem->Position.Angle(neigh->Position);
        }
        sum_angle += sum_angle_neigh / elem->Neighbours.size();
    }
    return sum_angle/detector.GetNChannels();
}

MCFakeReconstructed::MCFakeReconstructed() :
    cb(find_detector<decltype(cb)>()),
    cb_avg_angle(get_avg_opening_angle(*cb)),
    pid(find_detector<decltype(pid)>()),
    taps(find_detector<decltype(taps)>()),
    taps_avg_angle(get_avg_opening_angle(*taps)),
    tapsveto(find_detector<decltype(tapsveto)>())
{

}



void find_closest_ch(const vec3& pos, const Detector_t& detector,
                     unsigned& min_ch, double& min_angle) {
    // assume to have at least one channel...
    min_ch = 0;
    min_angle = pos.Angle(detector.GetPosition(min_ch));
    for(unsigned ch=1;ch<detector.GetNChannels();ch++) {
        const double angle = pos.Angle(detector.GetPosition(ch));
        if(angle < min_angle) {
            min_ch = ch;
            min_angle = angle;
        }
    }
}


bool do_calo_veto(const TParticle& p,
                  const ClusterDetector_t& calo,
                  double max_avg_angle,
                  const Detector_t& veto,
                  TEventData& data)
{
    unsigned calo_ch;
    double calo_angle;
    find_closest_ch(p.p, calo, calo_ch, calo_angle);

    // Calo determines if this particle could be detected
    if(calo_angle > max_avg_angle)
        return false;

    data.Clusters.emplace_back(
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
        unsigned veto_ch;
        double veto_angle;
        find_closest_ch(p.p, veto, veto_ch, veto_angle);
        data.Clusters.emplace_back(
                    veto.GetPosition(veto_ch).Unit(),
                    vetoE,
                    0, // timing zero for now
                    veto.Type,
                    veto_ch
                    );
        type |= veto.Type;
    }

    data.Candidates.emplace_back(
                type,
                p.Ek(),
                p.p.Theta(),
                p.p.Phi(),
                0,
                0, // cluster size 0 for now...
                vetoE,
                std_ext::NaN, // no tracker energy
                vetoE>0 ? TClusterList{std::prev(data.Clusters.end(), 2), std::prev(data.Clusters.end(), 1)}
                        : TClusterList{std::prev(data.Clusters.end())}
                );

    data.Particles.Add(std::make_shared<TParticle>(p.Type(), std::prev(data.Candidates.end())));

    return true;
}

const TEventData& MCFakeReconstructed::Get(const TEventData& mctrue)
{
    dataptr = std_ext::make_unique<TEventData>(mctrue.ID);
    TEventData& data = *dataptr;

    for(const TParticlePtr& p : mctrue.Particles.GetAll()) {
        // prefer CB to speed up matching
        if(!do_calo_veto(*p, *cb, cb_avg_angle, *pid, data))
            do_calo_veto(*p, *taps, taps_avg_angle, *tapsveto, data);
    }

    // data contains now "perfect" candidates/particles


    // some simple trigger stuff
    data.Trigger.CBTiming = 0;


    LOG(INFO) << data;

//    for(const auto& cand : data.Candidates) {
//        cout << *cand.FindCaloCluster() << endl;
//    }

    return data;
}