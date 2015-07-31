#include "Energy.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"
#include "calibration/CalibrationDataManager.h"

#include "tree/TDetectorRead.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

Energy::Energy(Detector_t::Type_t detectorType,
               std::shared_ptr<CalibrationDataManager> calmgr,
               Calibration::Converter::ptr_t converter,
               double defaultPedestal,
               double defaultGain,
               double defaultThreshold,
               double defaultRelativeGain) :
    Calibration::PhysicsModule(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_Energy"
           ),
    DetectorType(detectorType),
    calibrationManager(calmgr),
    Converter(move(converter)),
    DefaultPedestal(defaultPedestal),
    Pedestals(),
    DefaultGain(defaultGain),
    Gains(),
    DefaultThreshold(defaultThreshold),
    Thresholds(),
    DefaultRelativeGain(defaultRelativeGain)
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");
}

void Energy::ApplyTo(const readhits_t& hits, extrahits_t& extrahits)
{
    const auto& dethits = hits.get_item(DetectorType);

    // now calibrate the Energies (ignore any other kind of hits)
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->GetChannelType() != Channel_t::Type_t::Integral)
            continue;



        // Values might already be filled
        // (for example by previous calibration run, or A2Geant unpacker),
        // then we apply the threshold and the relative gain only
        std::vector<double> values(0);

        // prefer RawData if available
        if(!dethit->RawData.empty()) {
            // the Converter is smart enough to account for reference Energys!
            values = Converter->Convert(dethit->RawData);

            // for pedestal calibration, we insert extra hits here
            // containing the raw values
            extrahits.emplace_back(
                        LogicalChannel_t{
                            dethit->GetDetectorType(),
                            Channel_t::Type_t::Pedestal,
                            dethit->Channel
                        },
                        values
                        );

            // apply pedestal/gain/threshold to each of the values (might be multihit)
            for(double& value : values) {
                if(Pedestals.empty())
                    value -= DefaultPedestal;
                else
                    value -= Pedestals[dethit->Channel];

                if(Gains.empty())
                    value *= DefaultGain;
                else
                    value *= Gains[dethit->Channel];
            }

        }
        else {
            // maybe the values are already filled
            values = dethit->Values;
            dethit->Values.resize(0);
        }

        // always apply the threshold cut and the relative gains
        dethit->Values.reserve(values.size());

        for(double value : values) {
            if(RelativeGains.empty())
                value *= DefaultRelativeGain;
            else
                value *= RelativeGains[dethit->Channel];

            const double threshold = Thresholds.empty()
                                     ? DefaultThreshold
                                     : Thresholds[dethit->Channel];
            if(value<threshold)
                continue;

            // only add if it passes the threshold
            dethit->Values.push_back(value);
        }

    }
}

Energy::~Energy()
{

}

