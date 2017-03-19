#include "Energy.h"

#include "calibration/DataManager.h"

#include "tree/TCalibrationData.h"
#include "tree/TDetectorReadHit.h"

#include "base/Logger.h"
#include "base/std_ext/math.h"
#include "base/math_functions/Linear.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

Energy::Energy(const detector_ptr_t& det,
               const std::shared_ptr<DataManager>& calmgr,
               const Calibration::Converter::ptr_t& converter,
               defaults_t defaultPedestals,
               defaults_t defaultGains,
               defaults_t defaultThresholds_Raw,
               defaults_t defaultThresholds_MeV,
               defaults_t defaultRelativeGains,
               Channel_t::Type_t channelType) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(det->Type)
        << "_"
        << ( channelType == Channel_t::Type_t::IntegralShort ? "Short" : "" )
        << "Energy"
           ),
    DetectorType(det->Type),
    ChannelType(channelType),
    calibrationManager(calmgr),
    Converter(move(converter)),
    Pedestals(det, "Pedestals", defaultPedestals),
    Gains(det, "Gains", defaultGains, "ggIM"),
    Thresholds_Raw(det, "Thresholds_Raw", defaultThresholds_Raw),
    Thresholds_MeV(det, "Thresholds_MeV", defaultThresholds_MeV),
    RelativeGains(det, "RelativeGains", defaultRelativeGains, "ggIM")
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");
}

Energy::~Energy()
{
}

void Energy::ApplyTo(const readhits_t& hits)
{
    const auto& dethits = hits.get_item(DetectorType);

    // now calibrate the Energies (ignore any other kind of hits)
    for(TDetectorReadHit& dethit : dethits) {
        if(dethit.ChannelType != ChannelType)
            continue;

        // prefer building from RawData if available
        if(!dethit.RawData.empty()) {
            // clear previously read values (if any)
            dethit.Values.resize(0);

            // apply pedestal/gain to each of the values (might be multihit)
            for(const double& conv : Converter->Convert(dethit.RawData)) {
                TDetectorReadHit::Value_t value(conv);
                value.Calibrated -= Pedestals.Get(dethit.Channel);

                const double threshold = Thresholds_Raw.Get(dethit.Channel);
                if(value.Calibrated<threshold)
                    continue;

                // calibrate with absolute gain
                value.Calibrated *= Gains.Get(dethit.Channel);

                dethit.Values.emplace_back(move(value));
            }
        }

        // apply relative gain and threshold on MC
        {
            auto it_value = dethit.Values.begin();
            while(it_value != dethit.Values.end()) {
                it_value->Calibrated *= RelativeGains.Get(dethit.Channel);

                if(IsMC) {
                    const double threshold = Thresholds_MeV.Get(dethit.Channel);
                    // erase from Values if below threshold
                    if(it_value->Calibrated<threshold) {
                        it_value = dethit.Values.erase(it_value);
                        continue;
                    }
                }

                ++it_value;
            }
        }
    }
}



std::list<Updateable_traits::Loader_t> Energy::GetLoaders()
{

    std::list<Updateable_traits::Loader_t> loaders;

    for(auto calibration : AllCalibrations) {

        auto loader = [this, calibration]
                (const TID& currPoint, TID& nextChangePoint)
        {
            TCalibrationData cdata;
            if(calibrationManager->GetData(
                   GetName()+"_"+ calibration->Name,
                   currPoint, cdata, nextChangePoint))
            {
                auto& values = calibration->Values;
                for (const auto& val: cdata.Data) {
                    if(values.size()<val.Key+1)
                        values.resize(val.Key+1);
                    values[val.Key] = val.Value;
                }

                // call notify load if present
                if(calibration->NotifyLoad)
                    calibration->NotifyLoad(*calibration);
            }
            else {
                LOG_IF(!calibration->Values.empty(), WARNING)
                        << "No calibration data found for " << calibration->Name
                        << " at changepoint TID=" << currPoint << ", using default values";
                calibration->Values.resize(0);
            }
        };

        loaders.emplace_back(loader);
    }

    return loaders;
}

void Energy::UpdatedTIDFlags(const TID& id)
{
    IsMC = id.isSet(TID::Flags_t::MC);
}


