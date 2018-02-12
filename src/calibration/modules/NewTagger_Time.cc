#include "NewTagger_Time.h"

#include "base/Logger.h"

#include <cstdint>

using namespace std;
using namespace ant;
using namespace ant::calibration;

NewTagger_Time::NewTagger_Time(shared_ptr<expconfig::detector::Tagger> tagg,
                   shared_ptr<DataManager> calmgr,
                   std::map<expconfig::detector::Tagger::TDCSector_t, Calibration::Converter::ptr_t> converters,
                   double defaultOffset,
                   shared_ptr<gui::PeakingFitFunction> fitFunction,
                   const interval<double>& timeWindow
                     ) :
    Time(tagg,
         calmgr,
         nullptr, // do not set any converter by default
         defaultOffset,
         fitFunction,
         timeWindow
         )
{
    // set each converter individually depending on the TDC-sector of the new tagger
    for(unsigned ch=0;ch<tagg->GetNChannels();ch++) {
        /// \todo better check if the given converters actually contain the sector
        Converters[ch] = converters[tagg->GetTDCSector(ch)];
    }

}
