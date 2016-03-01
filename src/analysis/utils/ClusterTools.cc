#include "ClusterTools.h"

#include "tree/TCluster.h"
#include "expconfig/ExpConfig.h"
#include "base/std_ext/math.h"

#include <vector>
#include <algorithm>

#include "TVector3.h"

using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis::utils;
using namespace std;

double ClusterTools::LateralMoment(const TCluster &cluster)
{
    
    const auto& hits = cluster.Hits;

    if(cluster.Hits.size() < 3)
        return NaN;
    
        
    const auto& center = cluster.Position;

    const auto setup = ant::ExpConfig::Setup::GetLastFound();

    if(!setup) {
        throw std::runtime_error("No Setup found");
    }

    const auto det = setup->GetDetector(cluster.DetectorType);

    const auto r0 = 7.0; // TODO FIXME

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

        const TVector3 d = (center - det->GetPosition(h->Channel));
        const auto r = d.Mag();
        ret += h->Channel * sqr(r);

    }

    return ret / (ret + sqr(r0)*(hit0->Energy + hit1->Energy));
}
