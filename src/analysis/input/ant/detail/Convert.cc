#include "Convert.h"
#include <memory>
#include "data/Event.h"
#include "tree/TEvent.h"
#include "tree/TCluster.h"
#include "tree/TTrack.h"
#include "tree/TTaggerHit.h"
#include "tree/TDataRecord.h"
#include "analysis/data/Track.h"
#include "analysis/Detector.h"

using namespace ant;
using namespace std;

shared_ptr<ant::Event> input::Convert(const TEvent &event)
{
    shared_ptr<Event> antevent = make_shared<Event>();

    for(auto& track : event.Tracks) {
        antevent->Reconstructed().Tracks().emplace_back(
                    Convert(track)
                    );
    }

    return std::move(antevent);
}


std::shared_ptr<Track> input::Convert(const TTrack &track)
{
    auto anttrack = make_shared<Track>(
                track.Energy,
                track.Theta,
                track.Phi,
                track.Time,
                0,
                ant::detector_t::None,
                track.VetoEnergy,
                0,
                0
                );
    return anttrack;
}
