#include "Convert.h"
#include <memory>
#include "data/Event.h"
#include "tree/TEvent.h"
#include "tree/TCluster.h"
#include "tree/TTrack.h"
#include "tree/TTaggerHit.h"
#include "tree/TDataRecord.h"
#include "analysis/data/Track.h"
#include "analysis/data/TaggerHit.h"
#include "analysis/Detector.h"
#include <limits>

using namespace ant;
using namespace ant::input;
using namespace std;

template <class Cont1, class Cont2>
void Copy(const Cont1& from, Cont2& to) {
    to.reserve(from.size());
    for(auto& from_element : from) {
        to.emplace_back(ant::input::Convert(from_element));
    }
}

shared_ptr<ant::Event> input::Convert(const TEvent &event)
{
    shared_ptr<Event> antevent = make_shared<Event>();

    Copy(event.Tracks,     antevent->Reconstructed().Tracks());
    Copy(event.TaggerHits, antevent->Reconstructed().TaggerHits());

    return std::move(antevent);
}


std::shared_ptr<Track> input::Convert(const TTrack &track)
{
    /// @todo implement cluster size
    /// @todo add clusters to ant::Track
    auto anttrack = make_shared<Track>(
                        track.Energy,
                        track.Theta,
                        track.Phi,
                        track.Time,
                        0,
                        ant::detector_t::None,
                        track.VetoEnergy,
                        track.TrackerEnergy
                        );
    return anttrack;
}

std::shared_ptr<TaggerHit> input::Convert(const TTaggerHit& taggerhit)
{
    ///@todo implement something for tagger channel
    auto anttaggerhit = make_shared<TaggerHit>(
                            0,
                            taggerhit.PhotonEnergy,
                            taggerhit.Time
                            );
    return anttaggerhit;
}
