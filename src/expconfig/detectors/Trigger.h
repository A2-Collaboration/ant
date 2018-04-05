#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"
#include "reconstruct/Reconstruct_traits.h"

#include <bitset>
#include <stdexcept>

namespace ant {
namespace expconfig {
namespace detector {
struct Trigger :
        Detector_t,
        UnpackerAcquConfig
{

    Trigger() : Detector_t(Detector_t::Type_t::Trigger) {}

    virtual vec3 GetPosition(unsigned) const override {
        // when you ask the trigger detector for positions,
        // this is certainly a bug :)
        throw Exception("The trigger detector knows nothing about positions.");
    }
    virtual unsigned GetNChannels() const override {
        throw Exception("The trigger detector knows nothing about number of channels.");
    }
    virtual void SetElementFlags(unsigned, const ElementFlags_t&) override {
        throw Exception("The trigger detector cannot set flags on elements.");
    }
    virtual const ElementFlags_t& GetElementFlags(unsigned) const override {
        // do not throw here as innocent people (such as Reconstruct) may still ask
        static ElementFlags_t none;
        return none;
    }

    struct ReferenceTimingHitMapping_t {
        LogicalChannel_t LogicalChannel;
        unsigned AcquRawChannel;
        ReferenceTimingHitMapping_t(unsigned logicalChannel, unsigned acquRawChannel) :
            LogicalChannel{Detector_t::Type_t::Trigger, Channel_t::Type_t::Timing, logicalChannel},
            AcquRawChannel(acquRawChannel)
        {}
        operator LogicalChannel_t() const {
            return LogicalChannel;
        }
    };

    static const ReferenceTimingHitMapping_t Reference_CATCH_TaggerCrate;
    static const ReferenceTimingHitMapping_t Reference_V1190_TaggerTDC1;
    static const ReferenceTimingHitMapping_t Reference_V1190_TaggerTDC2;
    static const ReferenceTimingHitMapping_t Reference_V1190_TaggerTDC3_1;
    static const ReferenceTimingHitMapping_t Reference_V1190_TaggerTDC3_2;
    static const ReferenceTimingHitMapping_t Reference_CATCH_CBCrate;

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    // Define some scaler names of commonly used scalers
    struct ScalerName {
        static const std::string TotalLivetime;
        static const std::string Exptrigger_1MHz;
        static const std::string ExpTrigger;
        static const std::string L1Trigger;

        static const std::string EPTReferenceOR;

        static const std::string PairSpecGate;
        static const std::string TaggerReferenceOR;
        static const std::string PbGlass;
        static const std::string Paddle;
        static const std::string IonChamber;
        static const std::string FaradayCup;
        static const std::string Beampolmon_1MHz;
    };

    // base class knows nothing about scaler references
    virtual std::string GetScalerReference(const std::string& scalername) const =0;

};


struct Trigger_2014 :
        Trigger,
        ReconstructHook::DetectorReadHits
{

    static const ReferenceTimingHitMapping_t Reference_V1190_TAPSPbWO4;

    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    virtual std::string GetScalerReference(const std::string& scalername) const override;

    virtual void ApplyTo(const readhits_t& hits) override;

    Trigger_2014() :
        patterns(9), // VUPROMs give nine 16bit values as trigger patterns
        scaler_mapping(MakeScalerMapping())
    {}

    /// \todo once more about the old trigger system is known,
    /// those methods could be generalized and/or put into base class (possibly virtual)

    /// \todo since those patterns might be needed in Reconstruct, it might be better NOT to
    /// use those methods too much

    std::bitset<16> GetL1Pattern() const;
    std::bitset<16> GetL2Pattern() const;
    std::bitset<64> GetMultiplicityPattern() const;
    std::uint16_t   GetMultiplicityValue() const;
    std::bitset<16> GetHelicityPattern() const;
    std::bitset<16> GetTriggerFiredPattern() const;

protected:
    std::vector<std::bitset<16>> patterns;

    using scaler_mapping_t = std::map<std::string, unsigned>;
    static scaler_mapping_t MakeScalerMapping();
    const scaler_mapping_t scaler_mapping;


};

struct Trigger_2007 : Trigger {
    using Trigger::Trigger;

    virtual std::string GetScalerReference(const std::string&) const override;
};

}}} // namespace ant::expconfig::detector
