#include "ClusterTools.h"

#include "tree/TCluster.h"
#include "expconfig/ExpConfig.h"
#include "base/std_ext/math.h"

#include <vector>
#include <algorithm>

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis::utils;
using namespace std;

ClusterTools::ClusterTools()
{
    auto& setup = ExpConfig::Setup::Get();
    for(const auto& detector : setup.GetDetectors()) {
        auto clusterdetector = dynamic_pointer_cast<ClusterDetector_t, Detector_t>(detector);
        if(clusterdetector)
            cluster_detectors.emplace(clusterdetector->Type, clusterdetector);
    }
}


ClusterTools::det_t::det_t(const det_ptr_t& det) : Detector(det)
{
    for(unsigned ch=0;ch<Detector->GetNChannels();ch++) {
        // build average of distances to neighbours
        const auto& element = Detector->GetClusterElement(ch);
        double r_sum = 0;
        for(const auto& n : element->Neighbours) {
            r_sum += (Detector->GetPosition(n) - element->Position).R();
        }
        R0.emplace_back(r_sum/element->Neighbours.size());
    }
}

double ClusterTools::LateralMoment(const TCluster &cluster) const
{

    const auto& hits = cluster.Hits;

    if(cluster.Hits.size() < 3)
        return NaN;

    auto it_det = cluster_detectors.find(cluster.DetectorType);
    if(it_det == cluster_detectors.end())
        return NaN;
    const auto& det = it_det->second;


    const auto& center = cluster.Position;

    const auto r0 = det.R0[cluster.CentralElement];

    auto hit = hits.cbegin();

    // two highest energy hits
    auto hit0 = hit++;
    auto hit1 = hit++;

    double ret = 0.0;

    for(; hit!=hits.cend(); ++hit) {

        auto h = hit;

        if(h->Energy > hit0->Energy) {
            swap(h, hit0);
        } else if(h->Energy > hit1->Energy) {
            swap(h,hit1);
        }

        const auto d = center - det.Detector->GetPosition(h->Channel);
        const auto r = d.R();
        ret += h->Channel * sqr(r);

    }

    return ret / (ret + sqr(r0)*(hit0->Energy + hit1->Energy));
}

double ClusterTools::EffectiveRadius(const TCluster &cluster) const
{
    const auto& hits = cluster.Hits;

    if(cluster.Hits.size() < 3)
        return NaN;

    auto it_det = cluster_detectors.find(cluster.DetectorType);
    if(it_det == cluster_detectors.end())
        return NaN;
    const auto& det = it_det->second;


    const auto& center = cluster.Position;

    double effR = 0., e = 0.;

    for (const auto& hit : hits) {
        const auto r = radian_to_degree(center.Angle(det.Detector->GetPosition(hit.Channel)));
        effR += hit.Energy * sqr(r);
        e += hit.Energy;
    }

    effR /= e;

    return sqrt(effR);
}

