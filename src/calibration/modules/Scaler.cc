#include "Scaler.h"

#include "tree/TDetectorRead.h"

#include "base/std_ext.h"
#include "base/Logger.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;


Scaler::Scaler(Detector_t::Type_t detectorType, Calibration::Converter::ptr_t converter) :
    Calibration::BaseModule(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_Scaler"
           ),
    DetectorType(detectorType),
    Converter(move(converter))
{}

Scaler::~Scaler() {}

void Scaler::ApplyTo(const readhits_t& hits, extrahits_t&)
{
    // search for to be calibrated scalers
    const auto& dethits = hits.get_item(DetectorType);

    // now calibrate the scalers using the Converter
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->GetChannelType() != Channel_t::Type_t::Scaler)
            continue;
        dethit->Values = Converter->Convert(dethit->RawData);
    }
}
