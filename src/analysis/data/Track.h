#pragma once

#include "base/printable.h"
#include "base/types.h"
#include "analysis/Detector.h"

#include <ostream>
#include <memory>
#include <vector>

namespace ant {

/**
 * @brief The Track class
 * Representation of GoAT information, with emphasis on
 * physical particle information, not on detector information!
 */
class Track: public ant::printable_traits
{
private:
    mev_t clusterEnergy;
    radian_t theta;
    radian_t phi;
    ns_t time;
    clustersize_t clusterSize;
    detector_t detector;
    mev_t vetoEnergy;
    mev_t trackerEnergy;
public:
    Track(const mev_t& _clusterEnergy,
          const radian_t& _theta,
          const radian_t& _phi,
          const ns_t& _time,
          const clustersize_t& _clusterSize,
          const detector_t& _detector,
          const mev_t& _vetoEnergy,
          const mev_t& _trackerEnergy
          ) :
        clusterEnergy(_clusterEnergy),
        theta(_theta),
        phi(_phi),
        time(_time),
        clusterSize(_clusterSize),
        detector(_detector),
        vetoEnergy(_vetoEnergy),
        trackerEnergy(_trackerEnergy)
    {}


    mev_t ClusterEnergy() const { return clusterEnergy; }
    radian_t Theta() const { return theta; }
    radian_t Phi() const { return phi; }
    ns_t Time() const { return time; }
    clustersize_t ClusterSize() const { return clusterSize; }
    detector_t Detector() const { return detector; }
    mev_t VetoEnergy() const { return vetoEnergy; }
    mev_t TrackerEnergy() const { return trackerEnergy; }

    virtual std::ostream &Print(std::ostream &stream) const override;

};

using TrackPtr  = std::shared_ptr<Track>;
using TrackList = std::vector<TrackPtr>;

}
