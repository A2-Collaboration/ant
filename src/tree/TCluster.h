#pragma once

#include "base/root_printable.h"
#include "base/Detector_t.h"
#include "base/std_ext/math.h"

#include "TVector3.h"

#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <memory>


namespace ant {

struct TClusterHit;
using TClusterHitList = std::vector<TClusterHit>;

struct TCluster;
using TClusterPtr = std::shared_ptr<TCluster>;
using TClusterList = std::vector<TClusterPtr>;

struct TClusterHit : printable_traits
{
    std::uint32_t Channel;
    double Energy = std::numeric_limits<double>::quiet_NaN();
    double Time = std::numeric_limits<double>::quiet_NaN();

    struct Datum
    {
        Channel_t::Type_t Type;
        double Value;

        Datum(Channel_t::Type_t type, double value) :
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

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TClusterHit Ch=" << Channel << ", Energy=" << Energy << ", Time=" << Time;
        for(const auto& datum : Data) {
            s << ", " << Channel_t::ToString(datum.Type) << "=" << datum.Value;
        }
        return s;
    }

    TClusterHit() {}
    virtual ~TClusterHit() {}
};

struct TCluster : printable_traits
{
    double Energy;
    double Time;
    TVector3 Position;
    Detector_t::Type_t DetectorType;
    std::uint32_t CentralElement;
    std::uint32_t Flags;
    double ShortEnergy;

    TClusterHitList Hits;


    TCluster(
            const TVector3& pos,
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

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TCluster: " << Hits.size() << " hits @" << Position
          << ", Energy=" << Energy << " ShortEnergy=" << ShortEnergy
          << " CentralElement=" << CentralElement
          << " Detector=" << Detector_t::ToString(DetectorType) << std::endl;
        for(const auto& hit : Hits)
            s << hit << std::endl;
        return s;
    }

    // the cluster alorithm may set those flags
    enum class Flags_t : std::uint8_t {
        TouchesHole, Split, Unmatched
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

    TCluster() : Energy(), Time(),
        Position(),
        DetectorType(), CentralElement(), Flags(), ShortEnergy() {}
    virtual ~TCluster() {}
};

}

