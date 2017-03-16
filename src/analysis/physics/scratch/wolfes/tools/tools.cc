#include "tools.h"

using namespace std;
using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::physics;

#include <numeric>


double tools::getChargedClusterE(const TClusterList& clusters)
{
    auto accE = 0.0;
    for(const auto& cl : clusters)
    {
        if(cl.DetectorType == Detector_t::Type_t::PID ||
           cl.DetectorType == Detector_t::Type_t::TAPSVeto)
            accE += cl.Energy;
    }
    return accE;
}

double tools::getChargedCandidateE(const TCandidateList& cands)
{
    if (cands.size() == 0) return 0;
    return accumulate(cands.begin(),cands.end(),
                      cands[0].VetoEnergy,
                      [] (double acc, const TCandidate& ca)
                      {
                          return acc + ca.VetoEnergy;
                      });
}

