#pragma once

#include "base/printable.h"

#include "TVector3.h"

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ant {

/**
 * @brief The Detector_t struct is the minimal base class for all detectors
 */
struct Detector_t : printable_traits {

    // changing types here breaks the file format
    // the only operation allowed is to add detector types!
    enum class Type_t : std::uint8_t {
        Trigger,
        Tagger, TaggerMicro, EPT, Moeller, PairSpec,
        CB, PID, MWPC0, MWPC1,
        TAPS, TAPSVeto, Cherenkov
    };

    // Any_t represents a collection of detectors
    struct Any_t : printable_traits {

        static Any_t None() { return Any_t(0); }
        static Any_t MWPC() { return None() | Type_t::MWPC0 | Type_t::MWPC1; }
        static Any_t CB()   { return MWPC() | Type_t::PID | Type_t::CB; }
        static Any_t TAPS() { return None() | Type_t::TAPS | Type_t::TAPSVeto; }
        static Any_t Veto() { return None() | Type_t::PID | Type_t::TAPSVeto; }

        Any_t(const Type_t& type);

        bool operator==(const Any_t& other) const;
        Any_t operator&(const Any_t& other) const;
        Any_t operator|(const Any_t& other) const;
        Any_t operator^(const Any_t& other) const;
        explicit operator bool() const;
        operator std::string() const;

        virtual std::ostream& Print(std::ostream& stream) const override;
    private:
        Any_t(std::uint64_t bitfield_) : bitfield(bitfield_) {}
        std::uint64_t bitfield;

    }; // struct Any_t

    static const char* ToString(const Type_t& type);

    const Type_t Type;

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
    virtual unsigned GetNChannels() const = 0;

    class Exception : std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };

    virtual ~Detector_t() = default;
    virtual std::ostream& Print(std::ostream& stream) const override {
        return stream << "Detector_t " << ToString(Type);
    }
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
        throw Exception("You cannot ask a TaggerDetector_t for its position");
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

inline std::ostream& Detector_t::Any_t::Print(std::ostream& stream) const  {
    typename std::underlying_type<Type_t>::type i = 0;
    decltype(bitfield) temp  = 1;
    while(temp<=bitfield) {
        if(bitfield & temp) {
            stream << Detector_t::ToString(
                          static_cast<Detector_t::Type_t>(i)
                          )
                   << " ";
        }
        ++i;
        temp <<= 1;
    }
    return stream;
}

inline Detector_t::Any_t::Any_t(const Type_t& type)  :
    bitfield(1 << static_cast<typename std::underlying_type<Type_t>::type>(type))
{}

inline bool Detector_t::Any_t::operator==(const Any_t& other) const {
    return this->bitfield == other.bitfield;
}
inline Detector_t::Any_t Detector_t::Any_t::operator&(const Any_t& other) const {
    return this->bitfield & other.bitfield;
}
inline Detector_t::Any_t Detector_t::Any_t::operator|(const Any_t& other) const {
    return this->bitfield | other.bitfield;
}
inline Detector_t::Any_t Detector_t::Any_t::operator^(const Any_t& other) const {
    return this->bitfield ^ other.bitfield;
}
inline Detector_t::Any_t::operator bool() const {
    return bitfield;
}
inline Detector_t::Any_t::operator std::string() const {
    std::stringstream s;
    Print(s);
    return s.str();
}

inline bool operator!=(const Detector_t::Any_t& any1, const Detector_t::Any_t& any2) {
    return !(any1 == any2);
}
inline Detector_t::Any_t& operator&=(Detector_t::Any_t& any, const Detector_t::Any_t& other) {
    any = any & other;
    return any;
}
inline Detector_t::Any_t& operator|=(Detector_t::Any_t& any, const Detector_t::Any_t& other) {
    any = any | other;
    return any;
}
inline Detector_t::Any_t& operator^=(Detector_t::Any_t& any, const Detector_t::Any_t& other) {
    any = any ^ other;
    return any;
}

inline bool operator==(const Detector_t::Type_t& type, const Detector_t::Any_t& any) {
    return any == type;
}
inline bool operator!=(const Detector_t::Type_t& type, const Detector_t::Any_t& any) {
    return any != type;
}
inline Detector_t::Any_t operator|(const Detector_t::Type_t& type1,
                                   const Detector_t::Type_t& type2) {
    return static_cast<Detector_t::Any_t>(type1) | type2;
}

inline bool operator&(const Detector_t::Type_t& type, const Detector_t::Any_t& any) {
    return static_cast<bool>(any & type);
}
inline bool operator^(const Detector_t::Type_t& type, const Detector_t::Any_t& any) {
    return static_cast<bool>(any ^ type);
}



} // namespace ant



