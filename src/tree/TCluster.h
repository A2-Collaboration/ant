#pragma once

#include "tree/TDetectorReadHit.h"

#include "base/Detector_t.h"
#include "base/std_ext/math.h"
#include "base/std_ext/shared_ptr_container.h"

#include "base/vec/vec3.h"

#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <memory>


namespace ant {

struct TClusterHit;
using TClusterHitList = std::vector<TClusterHit>;

struct TCluster;
using TClusterList = std_ext::shared_ptr_container<TCluster>;
using TClusterPtr = std_ext::cc_shared_ptr<TCluster>;
using TClusterPtrList = std::vector<TClusterPtr>;

struct TClusterHit
{
    std::uint32_t Channel;
    double Energy = std_ext::NaN;
    double Time = std_ext::NaN;

    struct Datum
    {
        Channel_t::Type_t         Type;
        TDetectorReadHit::Value_t Value;

        Datum(Channel_t::Type_t type,
              const TDetectorReadHit::Value_t& value) :
            Type(type),
            Value(value)
        {}

        template<class Archive>
        void serialize(Archive& archive) {
            archive(Type, Value);
        }

        Datum() {}
    };

    std::vector<Datum> Data;

    TClusterHit(unsigned channel,
                double energy,
                double time) :
        Channel(channel),
        Energy(energy),
        Time(time)
    {
        static_assert(sizeof(Channel)>=sizeof(channel),
                      "Parameter channel does not fit into TClusterHit::Channel");
    }

    bool IsSane() const {
        return std::isfinite(Time) && std::isfinite(Energy);
    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Channel, Energy, Time, Data);
    }

    friend std::ostream& operator<<( std::ostream& s, const TClusterHit& o) {
        s << "TClusterHit Ch=" << o.Channel << ", Energy=" << o.Energy << ", Time=" << o.Time;
        for(const auto& datum : o.Data) {
            s << ", " << Channel_t::ToString(datum.Type) << "=" << datum.Value;
        }
        return s;
    }

    TClusterHit() = default;
};

struct TCluster
{
    double Energy;
    double Time;
    vec3 Position;
    Detector_t::Type_t DetectorType;
    std::uint32_t CentralElement;
    std::uint32_t Flags;
    double ShortEnergy;

    TClusterHitList Hits;


    TCluster(
            const vec3& pos,
            double E,
            double t,
            const Detector_t::Type_t& type,
            const unsigned central,
            const std::vector<TClusterHit>& hits = {}
            ) :
        Energy(E),
        Time(t),
        Position(pos),
        DetectorType(type),
        CentralElement(central),
        Flags(0),
        ShortEnergy(std_ext::NaN),
        Hits(hits)
    {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Energy, Time, Position,
                DetectorType, CentralElement, Flags, ShortEnergy, Hits);
    }

    friend std::ostream& operator<<( std::ostream& s, const TCluster& o) {
        s << "TCluster: " << o.Hits.size() << " hits @" << o.Position
          << ", Energy=" << o.Energy << " ShortEnergy=" << o.ShortEnergy
          << " CentralElement=" << o.CentralElement
          << " Detector=" << Detector_t::ToString(o.DetectorType) << std::endl;
        for(const auto& hit : o.Hits)
            s << hit << std::endl;
        return s;
    }

    /**
     * @brief The Flags_t enum set by the Reconstruction algorithm
     *
     * TouchesHoleCrystal: At least one crystal of cluster touches hole (ignored element)
     * TouchesHoleCentral: Central element (highest energy) touches hole
     * Split: The cluster was splitted from a larger number of hits (see clustering algorithm)
     * Unmatched: The cluster was not assigned to any candidate
     */
    enum class Flags_t : std::uint8_t {
        TouchesHoleCrystal, TouchesHoleCentral, Split, Unmatched
    };

    void SetFlag(Flags_t flag, bool set = true) {
        const std::uint32_t mask = 1 << static_cast<std::uint8_t>(flag);
        if(set) {
            Flags |= mask;
        }
        else {
            Flags &= ~mask;
        }
    }
    bool HasFlag(Flags_t flag) const {
        const std::uint32_t mask = 1 << static_cast<std::uint8_t>(flag);
        return (Flags & mask) != 0;
    }

    double GetPSARadius() const {
        using namespace ant::std_ext;
        return std::sqrt(sqr(Energy) + sqr(ShortEnergy));

    }
    double GetPSAAngle() const {
        return std::atan2(ShortEnergy, Energy);
    }

    bool isSane() const {
        return std::isfinite(Energy) && std::isfinite(Time);
    }

    TCluster() : Energy(std_ext::NaN), Time(std_ext::NaN),
        Position(),
        DetectorType(), CentralElement(), Flags(), ShortEnergy()
    {}

    TCluster(const TCluster&) = delete;
    TCluster(TCluster&&) = delete;
    TCluster& operator=(const TCluster&) = delete;
    TCluster& operator=(TCluster&&) = delete;

};

}

