#pragma once

#include "base/printable.h"

#include "base/vec/vec3.h"
#include "base/std_ext/math.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <string>

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
        TAPS, TAPSVeto, Cherenkov,
        Raw
    };

    // Any_t represents a collection of detectors
    struct Any_t : printable_traits {
        static const Any_t None;
        static const Any_t Tracker; // i.e. MWPC
        static const Any_t CB_Apparatus; // i.e. PID, MWPC, CB
        static const Any_t TAPS_Apparatus; // i.e. TAPS, TAPSVeto
        static const Any_t Calo; // i.e. CB or TAPS calorimeter
        static const Any_t Veto; // i.e. PID or TAPSVeto

        Any_t(const Type_t& type);

        bool operator==(const Any_t& other) const;
        Any_t operator&(const Any_t& other) const;
        Any_t operator|(const Any_t& other) const;
        Any_t operator^(const Any_t& other) const;
        explicit operator bool() const;
        operator std::string() const;

        virtual std::ostream& Print(std::ostream& stream) const override;

        // support cereal
        template<class Archive>
        void serialize(Archive archive) {
            archive(bitfield);
        }

    private:
        Any_t(std::uint64_t bitfield_) : bitfield(bitfield_) {}
        std::uint64_t bitfield;

    }; // struct Any_t

    static const char* ToString(const Type_t& type);
    static Type_t      FromString(const std::string& str);

    const Type_t Type;

    // Element_t is the minimum information,
    // derived classes may extend this
    struct Element_t {
        Element_t(unsigned channel, const vec3& position) :
            Channel(channel),
            Position(position)
        {}
        unsigned Channel; // unique within Detector for all time!
        vec3 Position;
    };

    virtual unsigned GetNChannels() const = 0;
    virtual vec3 GetPosition(unsigned channel) const = 0;
    virtual void SetIgnored(unsigned channel) = 0;
    virtual bool IsIgnored(unsigned channel) const = 0;

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
        BitPattern, Raw
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
                const vec3& position,
                const std::vector<unsigned>& neighbours,
                double moliereRadius,
                double criticalE,
                double radiationLength,
                bool touchesHole = false
                ) :
            Detector_t::Element_t(channel, position),
            Neighbours(neighbours),
            MoliereRadius(moliereRadius),
            CriticalE(criticalE),
            RadiationLength(radiationLength),
            TouchesHole(touchesHole)
        {}
        std::vector<unsigned> Neighbours;
        double MoliereRadius;
        double CriticalE;       // in MeV
        double RadiationLength; // X0 in cm
        bool TouchesHole;
    };

    virtual const Element_t* GetClusterElement(unsigned channel) const = 0;

protected:
    ClusterDetector_t(const Type_t& type) :
        Detector_t(type) {}
};

struct TaggerDetector_t : Detector_t {

    virtual double GetPhotonEnergy(unsigned channel) const = 0;

    /**
     * @brief GetPhotonEnergyWidth
     * @param channel
     * @return Width of channel in MeV
     *
     * Calculates the with from the distance to the neighbor channels.
     * Override if this is not what you want for your tagger.
     *
     * @note Only works for taggers with more than one channel.
     */
    virtual double GetPhotonEnergyWidth(unsigned channel) const;
    bool TryGetChannelFromPhoton(double photonEnergy, unsigned& channel) const;

    struct taggeff_t
    {
        double Value;
        double Error;
        taggeff_t(double value = std_ext::NaN, double error = std_ext::NaN):
            Value(value),
            Error(error){}
    };
    virtual taggeff_t GetTaggEff(unsigned channel) const = 0;
    virtual void SetTaggEff(unsigned channel, const taggeff_t& taggEff) = 0;


    virtual vec3 GetPosition(unsigned) const final {
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
        unsigned  Channel;
        double    ElectronEnergy;
        taggeff_t TaggEff;
    };

    double BeamEnergy;

    TaggerDetector_t(const Type_t& type,
                     double beamEnergy
                     ) :
        Detector_t(type),
        BeamEnergy(beamEnergy)
    {}
};

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



