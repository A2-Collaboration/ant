#pragma once


namespace ant {

class TCluster;

namespace analysis {
namespace utils {

struct ClusterTools {

    /**
     * @brief Calculate the Lateral Moment of a cluster
     * @param cluster
     * @return
     *
     * \f[
     *      lat = \frac{\sum_{i=2}^{n} E_i r_i^2 }{ (\sum_{i=2}^{n} E_i r_i^2 + r_0^0 E_0 + r_0^2 E_1) }
     * \f]
     * where \f$ E_{0} \f$ and \f$ E_{1} \f$ are the two hits highest energies.
     *
     * Reimplemented after BaBar
     */
     static double LateralMoment(const ant::TCluster& cluster);
};

}}}
