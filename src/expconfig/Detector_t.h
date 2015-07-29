#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <type_traits>

#include "TVector3.h"

namespace ant {

/**
 * @brief The Detector_t struct is the minimal base class for all detectors
 */
struct Detector_t {

    // changing types here breaks the file format
    // the only operation allowed is to add detector types!
    enum class Type_t : std::uint8_t {
        Trigger,
        Tagger, TaggerMicro, EPT, Moeller, PairSpec,
        CB, PID, MWPC0, MWPC1,
        TAPS, TAPSVeto, Cherenkov
    };
    using basetype = typename std::underlying_type<Type_t>::type;

    const Type_t Type;
    static const char* ToString(const Type_t& type);

    static constexpr std::uint32_t ToBitfield(const Type_t& type) {
        return 1 << static_cast<basetype>(type);
    }


    // Element_t is the minimum information,
    // derived classes may extend this
    struct Element_t {
        Element_t(unsigned channel, const TVector3& position) :
            Channel(channel),
            Position(position)
        {}
        unsigned Channel; // unique within Detector for all time!
        TVector3 Position;
    };

    virtual TVector3 GetPosition(unsigned channel) const = 0;

    virtual ~Detector_t() = default;

protected:
    Detector_t(const Type_t& type) :
        Type(type) {}
    Detector_t(const Detector_t&) = delete; // disable copy
};

/**
 * @brief The Channel_t struct enumerates the different channel types
 */
struct Channel_t {
    enum class Type_t : std::uint8_t {
        Timing,
        Integral, IntegralShort,
        IntegralAlternate, IntegralShortAlternate,
        BitPattern, Scaler, Pedestal, Counter
    };
    static bool IsIntegral(const Type_t& t);
    static const char* ToString(const Type_t& type);
};

/**
 * @brief The LogicalChannel_t struct uniquely identifies detector element within setup
 */
struct LogicalChannel_t {
    Detector_t::Type_t DetectorType;
    Channel_t::Type_t ChannelType;
    unsigned Channel;
};

struct ClusterDetector_t : Detector_t {

    struct Element_t : Detector_t::Element_t {
        Element_t(
                unsigned channel,
                const TVector3& position,
                const std::vector<unsigned>& neighbours,
                double moliereRadius
                ) :
            Detector_t::Element_t(channel, position),
            Neighbours(neighbours),
            MoliereRadius(moliereRadius)
        {}
        std::vector<unsigned> Neighbours;
        double MoliereRadius;
    };

    virtual const Element_t* GetClusterElement(unsigned channel) const = 0;

protected:
    ClusterDetector_t(const Type_t& type) :
        Detector_t(type) {}
};

struct TaggerDetector_t : Detector_t {

    virtual double GetPhotonEnergy(unsigned channel) const = 0;

    virtual bool TryGetChannelFromPhoton(double photonEnergy, unsigned& channel) const = 0;

    virtual TVector3 GetPosition(unsigned) const final {
        throw std::runtime_error("You cannot ask a TaggerDetector_t for its position");
    }


protected:
    // Tagger elements don't derive from Detector_t::Element
    // because positions are not meaningful for them (at least now)
    struct Element_t {
        Element_t(unsigned channel, double electronEnergy) :
            Channel(channel),
            ElectronEnergy(electronEnergy)
        {}
        unsigned Channel;
        double   ElectronEnergy;
    };

    double BeamEnergy;
    double ElectronEnergyWidth;

    TaggerDetector_t(const Type_t& type,
                     double beamEnergy,
                     double electronEnergyWidth
                     ) :
        Detector_t(type),
        BeamEnergy(beamEnergy),
        ElectronEnergyWidth(electronEnergyWidth)
    {}
};


inline const char* Channel_t::ToString(const Type_t& type)
{
    switch(type) {
    case Channel_t::Type_t::BitPattern:
        return "BitPattern";
    case Channel_t::Type_t::Counter:
        return "Counter";
    case Channel_t::Type_t::Integral:
        return "Integral";
    case Channel_t::Type_t::IntegralAlternate:
        return "IntegralAlternate";
    case Channel_t::Type_t::IntegralShort:
        return "IntegralShort";
    case Channel_t::Type_t::IntegralShortAlternate:
        return "IntegralShortAlternate";
    case Channel_t::Type_t::Scaler:
        return "Scaler";
    case Channel_t::Type_t::Timing:
        return "Timing";
    case Channel_t::Type_t::Pedestal:
        return "Pedestal";
    }
    throw std::runtime_error("Not implemented");
}

inline const char* Detector_t::ToString(const Type_t &type)
{
    switch(type) {
    case Detector_t::Type_t::CB :
        return "CB";
    case Detector_t::Type_t::Cherenkov:
        return "Cherenkov";
    case Detector_t::Type_t::MWPC0:
        return "MWPC0";
    case Detector_t::Type_t::MWPC1:
        return "MWPC1";
    case Detector_t::Type_t::PID:
        return "PID";
    case Detector_t::Type_t::Tagger:
        return "Tagger";
    case Detector_t::Type_t::TaggerMicro:
        return "TaggerMicro";
    case Detector_t::Type_t::EPT:
        return "EPT";
    case Detector_t::Type_t::Moeller:
        return "Moeller";
    case Detector_t::Type_t::PairSpec:
        return "PairSpec";
    case Detector_t::Type_t::TAPS:
        return "TAPS";
    case Detector_t::Type_t::TAPSVeto:
        return "TAPSVeto";
    case Detector_t::Type_t::Trigger:
        return "Trigger";
    }
    throw std::runtime_error("Not implemented");
}

inline bool Channel_t::IsIntegral(const Channel_t::Type_t& t) {
    switch(t) {
    case Type_t::Integral:
    case Type_t::IntegralShort:
    case Type_t::IntegralAlternate:
    case Type_t::IntegralShortAlternate:
        return true;
    default:
        return false;
    }
}

} // namespace ant
