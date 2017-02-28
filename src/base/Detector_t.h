#pragma once

#include "base/printable.h"
#include "base/bitflag.h"
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
    struct Any_t : bitflag<Type_t>, printable_traits {

        constexpr Any_t(Type_t t) : bitflag<Type_t>(t) {}
        constexpr Any_t(const bitflag<Type_t>& t) : bitflag<Type_t>(t) {}

        static const Any_t None;
        static const Any_t Tracker; // i.e. MWPC
        static const Any_t CB_Apparatus; // i.e. PID, MWPC, CB
        static const Any_t TAPS_Apparatus; // i.e. TAPS, TAPSVeto
        static const Any_t Calo; // i.e. CB or TAPS calorimeter
        static const Any_t Veto; // i.e. PID or TAPSVeto

        virtual std::ostream& Print(std::ostream& stream) const override;
    private:
        constexpr Any_t() = default;

    }; // struct Any_t

    // define the type of the detector, unique identifier!
    static const char* ToString(const Type_t& type);
    static Type_t      FromString(const std::string& str);
    const Type_t Type;


    enum class ElementFlag_t {
        Broken, // there's really nothing to do
        BadTDC, // no timing but energy could be used
        NoCalib // cannot be calibrated
    };

    // Element_t is the minimum information,
    // derived classes may (and will) extend this
    struct Element_t {
        Element_t(unsigned channel, const vec3& position) :
            Channel(channel),
            Position(position)
        {}
        unsigned Channel; // unique within Detector for all time!
        vec3 Position;
        bitflag<ElementFlag_t> Flags;
    };

    virtual unsigned GetNChannels() const = 0;
    virtual vec3 GetPosition(unsigned channel) const = 0;

    virtual void SetElementFlag(unsigned channel, ElementFlag_t flag) = 0;
    virtual void SetElementFlag(const std::vector<unsigned>& channels, ElementFlag_t flag);
    virtual bool HasElementFlag(unsigned channel, ElementFlag_t flag) const = 0;

    // common exception class
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
        bitflag<ElementFlag_t> Flags;
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


} // namespace ant



