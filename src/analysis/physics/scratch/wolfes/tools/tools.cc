#include "tools.h"

#include "base/std_ext/math.h"

#include <numeric>

using namespace std;
using namespace ant;
using namespace ant::std_ext;
using namespace ant::analysis;
using namespace ant::analysis::physics;

const double strictVetoCutAngle  = 15.0;

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

double tools::getCorrVetoEnergy(const TCandidate& photon, const TCandidate& proton)
{
    if (fabs(vec2::Phi_mpi_pi(proton.Phi - photon.Phi - M_PI)) > degree_to_radian(strictVetoCutAngle))
        return 0;
    return photon.VetoEnergy;
}


double tools::getPhotonVetoEnergy(const tools::protonSelection_t& sel, const bool strict)
{
    double acc = 0;
    for (const auto& ph: sel.Photons)
    {
        acc += strict ? getCorrVetoEnergy(*ph->Candidate,*sel.Proton->Candidate) : ph->Candidate->VetoEnergy;
    }
    return acc;
}

