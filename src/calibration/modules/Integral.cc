#include "Integral.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDetectorRead.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

Integral::Integral(Detector_t::Type_t detectorType,
        Calibration::Converter::ptr_t converter, const double defaultPedestal,
        const double defaultGain,
        const double defaultThreshold
        ) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_Integral"
           ),
    DetectorType(detectorType),
    Converter(move(converter)),
    DefaultPedestal(defaultPedestal),
    Pedestals(),
    DefaultGain(defaultGain),
    Gains(),
    DefaultThreshold(defaultThreshold),
    Thresholds()
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");
}

void Integral::ApplyTo(const readhits_t& hits, extrahits_t& extrahits)
{
    // search for to be calibrated Integrals
    const auto it_dethits = hits.find(DetectorType);
    if(it_dethits == hits.end())
        return;

    const auto& dethits = it_dethits->second;

    // now calibrate the Integrals (ignore any other kind of hits)
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->GetChannelType() != Channel_t::Type_t::Integral)
            continue;
        // the Converter is smart enough to account for reference Integrals!
        const std::vector<double>& values = Converter->Convert(dethit->RawData);
        dethit->Values.reserve(values.size());

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
        for(double value : values) {
            if(Pedestals.empty())
                value -= DefaultPedestal;
            else
                value -= Pedestals[dethit->Channel];

            if(Gains.empty())
                value *= DefaultGain;
            else
                value *= Gains[dethit->Channel];

            const double threshold = Thresholds.empty()
                    ? DefaultThreshold : Thresholds[dethit->Channel];
            if(value<threshold)
                continue;

            // only add if it passes the threshold
            dethit->Values.push_back(value);
        }
    }
}

