#include "TimeSmearingHack.h"

#include "TRandom.h"
#include "base/Logger.h"

using namespace ant;
using namespace ant::analysis;


double utils::TimeSmearingHack::CalculateTimeDifference(const TCandidatePtr& p1, const TCandidatePtr& p2) const
{
    auto time_diff = p1->Time - p2->Time;

    if(smearing_enabled) {
        if(   (p1->Detector & Detector_t::Any_t::CB_Apparatus   && p2->Detector & Detector_t::Any_t::TAPS_Apparatus)
              || (p1->Detector & Detector_t::Any_t::TAPS_Apparatus && p2->Detector & Detector_t::Any_t::CB_Apparatus  )) {

            return gRandom->Gaus(time_diff, cb_taps_sigma);
        }

        if(   p1->Detector & Detector_t::Any_t::CB_Apparatus
              && p2->Detector & Detector_t::Any_t::CB_Apparatus) {
            return gRandom->Gaus(time_diff, cb_cb_sigma);
        }

        if(   p1->Detector & Detector_t::Any_t::TAPS_Apparatus
              && p2->Detector & Detector_t::Any_t::TAPS_Apparatus) {
            return gRandom->Gaus(time_diff, taps_taps_sigma);
        }
    }

    return time_diff;

}

double utils::TimeSmearingHack::CalculateTimeDifference(const TCandidatePtr& p, const TTaggerHit& taggerhit) const
{
    auto time_diff = p->Time - taggerhit.Time;

    if(smearing_enabled) {
        if(p->Detector & Detector_t::Any_t::CB_Apparatus) {
            return gRandom->Gaus(time_diff, tagg_cb_sigma);
        }

        if(p->Detector & Detector_t::Any_t::TAPS_Apparatus) {
            return gRandom->Gaus(time_diff, tagg_taps_sigma);
        }
    }

    return time_diff;

}

double utils::TimeSmearingHack::GetTime(const TCandidatePtr& p) const
{
    const auto time = p->Time;

    if(p->Detector & Detector_t::Any_t::CB_Apparatus) {
        return gRandom->Gaus(time, cb_sigma);
    }

    if(p->Detector & Detector_t::Any_t::TAPS_Apparatus) {
        return gRandom->Gaus(time, taps_sigma);
    }

    return time;
}

utils::TimeSmearingHack::TimeSmearingHack() {
    LOG(WARNING) << "Ugly Time Smearing Hack in use!";
}
