#pragma once

#include "base/Detector_t.h"
#include "unpacker/UnpackerAcqu.h"
#include "reconstruct/Reconstruct_traits.h"

#include <stdexcept>

namespace ant {
namespace expconfig {
namespace detector {
struct Trigger :
        Detector_t,
        UnpackerAcquConfig,
        ReconstructHook::EventData
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
    virtual void SetIgnored(unsigned) override {
        throw Exception("The trigger detector cannot ignore channels.");
    }
    virtual bool IsIgnored(unsigned) const override {
        throw Exception("The trigger detector knows nothing about ignored channels.");
    }



    const LogicalChannel_t Reference_CATCH_TaggerCrate = {Type, Channel_t::Type_t::Timing, 1000};
    const LogicalChannel_t Reference_CATCH_CBCrate = {Type, Channel_t::Type_t::Timing, 1001};

    // for UnpackerAcquConfig
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    // for ReconstructHook::EventData
    // calculates the CBESum
    virtual void ApplyTo(TEventData& reconstructed) override;

    // Define some scaler names of commonly used scalers
    struct ScalerName {
        static const std::string Exptrigger_1MHz;
        static const std::string Beampolmon_1MHz;
        static const std::string TotalLivetime;
        static const std::string FaradayCup;
        static const std::string IonChamber;
        static const std::string PbGlass;
    };

    // base class knows nothing about scaler references
    virtual std::string GetScalerReference(const std::string& scalername) const =0;

};


struct Trigger_2014 : Trigger {


    const LogicalChannel_t Reference_V1190_TAPSPbWO4 = {Type, Channel_t::Type_t::Timing, 1002};



    virtual bool Matches(const TID& tid) const override;
    virtual void BuildMappings(
            std::vector<hit_mapping_t>&,
            std::vector<scaler_mapping_t>&) const override;

    virtual std::string GetScalerReference(const std::string& scalername) const override;

    Trigger_2014() : scaler_mapping(MakeScalerMapping()) {}

protected:
    using scaler_mapping_t = std::map<std::string, unsigned>;
    static scaler_mapping_t MakeScalerMapping();
    const scaler_mapping_t scaler_mapping;
};

}}} // namespace ant::expconfig::detector
