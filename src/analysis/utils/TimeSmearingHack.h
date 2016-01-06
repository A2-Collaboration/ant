#pragma once

#include "data/Particle.h"
#include "data/TaggerHit.h"

namespace ant {
namespace analysis {
namespace utils {

struct TimeSmearingHack {

    //extraced from data 2014-10-EPT_Prod
    double cb_cb_sigma     = 2.3;
    double taps_taps_sigma = 0.4;
    double cb_taps_sigma   = 2.5;
    double tagg_cb_sigma   = 1.55;
    double tagg_taps_sigma = 0.52;

    double cb_sigma        = 2.54;
    double taps_sigma      = 2.7;

    bool smearing_enabled = true;

    /**
     * @brief Calculate Time Difference of two candidates, and smear it accordingly if smearing is enabled
     * @param p1
     * @param p2
     * @return p1 time - p2 time
     */
    double CalculateTimeDifference(const data::CandidatePtr& p1, const data::CandidatePtr& p2) const;

    /**
     * @brief Calculate Time Difference of a candidate and a tagger hit, and smear it accordingly if smearing is enabled
     * @param p
     * @param taggerhit
     * @return time p - time taggerhit
     */
    double CalculateTimeDifference(const data::CandidatePtr& p , const data::TaggerHit& taggerhit) const;

    /**
     * @brief Get the Time inforamtion from a candidate, smeared if enabeld
     * @param p
     * @return time of p
     */
    double GetTime(const data::CandidatePtr& p) const;

    TimeSmearingHack();

};
}
}
}
