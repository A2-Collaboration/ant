#include "Time.h"
#include "analysis/plot/HistogramFactories.h"
#include "analysis/data/Event.h"
#include "analysis/utils/combinatorics.h"

#include "tree/TDetectorRead.h"

#include "base/Logger.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::analysis;
using namespace ant::analysis::data;

Time::Time(Detector_t::Type_t detectorType,
               Calibration::Converter::ptr_t converter,
        double defaultOffset,
        const interval<double>& timeWindow, // default {-inf, inf}
        const double defaultGain, // default gain is 1.0
        const std::vector< TKeyValue<double> >& gains
        ) :
    Calibration::Module(
        std_ext::formatter()
        << Detector_t::ToString(detectorType)
        << "_Time"
           ),
    DetectorType(detectorType),
    Converter(move(converter)),
    TimeWindow(timeWindow),
    DefaultOffset(defaultOffset),
    Offsets(),
    DefaultGain(defaultGain),
    Gains()
{
    if(Converter==nullptr)
        throw std::runtime_error("Given converter should not be nullptr");

    // fill a gain vector from given key-value pairs
    // for faster access (if some are given at all)
    if(gains.empty())
        return;
    unsigned maxkey = 0;
    for(const auto& gain : gains)
        maxkey = gain.Key>maxkey ? gain.Key : maxkey;
    Gains.resize(maxkey+1, DefaultGain);
    for(const auto& gain : gains)
        Gains[gain.Key] = gain.Value;
}

std::vector<std::list<TID> > Time::GetChangePoints() const {
   // vector<list<TID>> changePointLists;

//    for (const auto& calib_name: { Pedestals.Name, Gains.Name, Thresholds.Name, RelativeGains.Name}) {
 //       changePointLists.push_back(calibrationManager->GetChangePoints(
  //                                     GUI_CalibType::ConstructName(GetName(), calib_name)
   //                                    )
    //                               );
  //  }
    return {};
}

void Time::Update(size_t index, const TID& id) {
    LOG(INFO) << GetName() << ": Update called for index " << index << " with TID=" << id;
}

void Time::ApplyTo(const readhits_t& hits, extrahits_t&)
{
    const auto& dethits = hits.get_item(DetectorType);

    // now calibrate the Times (ignore any other kind of hits)
    for(TDetectorReadHit* dethit : dethits) {
        if(dethit->GetChannelType() != Channel_t::Type_t::Timing)
            continue;

        // the Converter is smart enough to account for reference Times!
        const auto& values = Converter->Convert(dethit->RawData);
        dethit->Values.reserve(values.size());

        // apply gain/offset to each of the values (might be multihit)
        for(double value : values) {
            if(Gains.empty())
                value *= DefaultGain;
            else
                value *= Gains[dethit->Channel];

            if(Offsets.empty())
                value -= DefaultOffset;
            else
                value -= Offsets[dethit->Channel];

            if(!TimeWindow.Contains(value))
                continue;

            dethit->Values.push_back(value);
        }
    }
}

Time::ThePhysics::ThePhysics(const string& name, const string& histName):
//                             unsigned nChannels0):
    Physics(name)
{
    hTime = HistFac.makeTH1D("Detectorname -Time",
                             "time [ns]","#",
                             BinSettings(200,-40,40),
                             histName);

}


void ant::calibration::Time::ThePhysics::ProcessEvent(const Event& event)
{
}

void ant::calibration::Time::ThePhysics::Finish()
{
}

void ant::calibration::Time::ThePhysics::ShowResult()
{
}
