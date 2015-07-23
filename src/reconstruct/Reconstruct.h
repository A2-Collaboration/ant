#ifndef ANT_RECONSTRUCT_H
#define ANT_RECONSTRUCT_H

#include <memory>
#include <list>

#include "Reconstruct_traits.h"
#include "tree/TEvent.h"



namespace ant {

class THeaderInfo;
class TDetectorRead;

namespace reconstruct {
class TrackBuilder;
}

class Reconstruct {
public:
    // You can only create the reconstruct machinery
    // if it's able to find its config. For this, it needs the
    // some THeaderInfo object
    Reconstruct(const THeaderInfo& headerInfo);

    // this method converts a TDetectorRead
    // into a calibrated TEvent
    std::unique_ptr<TEvent> DoReconstruct(TDetectorRead& detectorRead);

    ~Reconstruct();

private:

    /**
     * @brief The HitWithEnergy_t struct
     * Stores a TClusterHit together with its energy, if
     * exactly one TDetectorReadHit was of channel type integral
     */
    struct HitWithEnergy_t {
        TClusterHit Hit;
        double Energy;
        void MaybeSetEnergy(const TDetectorReadHit* readhit);
        HitWithEnergy_t(const TDetectorReadHit* readhit,
                        const std::vector<TClusterHitDatum>&& data);
    };

    template<typename T>
    using sorted_bydetectortype_t = std::map<Detector_t::Type_t, std::list< T > >;

    void BuildHits(
            sorted_bydetectortype_t<TDetectorReadHit*>&& sorted_readhits,
            sorted_bydetectortype_t<HitWithEnergy_t>& sorted_clusterhits,
            TTagger& event_tagger
            );
    void BuildClusters(

            );

    template<typename T>
    using shared_ptr_list = std::list< std::shared_ptr<T> >;

    using sorted_detectors_t = std::map<Detector_t::Type_t, std::shared_ptr<Detector_t> >;


    shared_ptr_list<CalibrationApply_traits>  calibrations;
    shared_ptr_list<Updateable_traits>        updateables;
    sorted_detectors_t                        sorted_detectors;

    std::unique_ptr<reconstruct::TrackBuilder> trackbuilder;

};

}

#endif // ANT_RECONSTRUCT_H
