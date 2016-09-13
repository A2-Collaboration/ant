#include "Tagger_QDC.h"

#include "tree/TDetectorReadHit.h"

using namespace std;
using namespace ant;
using namespace ant::calibration;

Tagger_QDC::Tagger_QDC(Detector_t::Type_t detectorType,
                       Calibration::Converter::ptr_t converter) :
    DetectorType(detectorType),
    Converter(converter)
{

}

Tagger_QDC::~Tagger_QDC()
{

}

void Tagger_QDC::ApplyTo(const ant::ReconstructHook::Base::readhits_t& hits)
{
    const auto& dethits = hits.get_item(DetectorType);
    for(TDetectorReadHit& dethit : dethits) {
        if(dethit.ChannelType != Channel_t::Type_t::Integral)
            continue;
        dethit.Converted = Converter->Convert(dethit.RawData);
        dethit.Values = dethit.Converted;
    }
}