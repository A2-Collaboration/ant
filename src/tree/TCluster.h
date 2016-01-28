#pragma once

#include "base/root_printable.h"
#include "base/Detector_t.h"
#include "base/std_ext/math.h"

#include "TVector3.h"

#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>


namespace ant {

struct TClusterHitDatum
{
    std::uint8_t Type;
    double Value;

    TClusterHitDatum(Channel_t::Type_t type, double value) :
        Type(static_cast<std::uint8_t>(type)),
        Value(value)
    {}

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Type, Value);
    }

    Channel_t::Type_t GetType() const {
        return static_cast<Channel_t::Type_t>(Type);
    }

    TClusterHitDatum() {}
};

struct TClusterHit : printable_traits
{
    std::uint32_t Channel;
    std::vector<TClusterHitDatum> Data;


    TClusterHit(unsigned channel,
                const std::vector<TClusterHitDatum>& data):
        Channel(channel),
        Data(data) {
        static_assert(sizeof(Channel)>=sizeof(channel),
                      "Parameter channel does not fit into TClusterHit::Channel");

    }

    template<class Archive>
    void serialize(Archive& archive) {
        archive(Channel, Data);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        s << "TClusterHit Ch=" << Channel << ": ";
        for(auto& datum : Data) {
            s << Channel_t::ToString(datum.GetType()) << "=" << datum.Value << " ";
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
    std::uint8_t DetectorType;
    std::uint32_t CentralElement;
    std::uint32_t Flags;
    double ShortEnergy;

    std::vector<TClusterHit> Hits;


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
        DetectorType(static_cast<std::uint8_t>(type)),
        CentralElement(central),
        Flags(0),
        ShortEnergy(std_ext::NaN),
        Hits(hits)
    {}

    template<class Archive>
    void load(Archive& archive) {
        double x,y,z;
        archive(Energy, Time,
                x, y, z,
                DetectorType, CentralElement, Flags, ShortEnergy, Hits);
        Position.SetXYZ(x,y,z);
    }

    template<class Archive>
    void save(Archive& archive) const {
        archive(Energy, Time,
                Position.X(), Position.Y(), Position.Z(),
                DetectorType, CentralElement, Flags, ShortEnergy, Hits);
    }

    Detector_t::Type_t GetDetectorType() const {
        return static_cast<Detector_t::Type_t>(DetectorType);
    }

    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TCluster: " << Hits.size() << " hits @" << Position
                 << ", Energy=" << Energy << " ShortEnergy=" << ShortEnergy
                 << " Central Element=" << CentralElement
                 << " Detector=" << Detector_t::ToString(GetDetectorType());
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

