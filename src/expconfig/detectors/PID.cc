#include "PID.h"
#include <cassert>
#include "base/std_ext.h"


#include "detail/PID_2004_elements.h"
#include "detail/PID_2009_05_elements.h"
#include "detail/PID_2009_06_elements.h"
#include "detail/PID_2009_07_elements.h"
#include "detail/PID_2014_elements.h"


using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;



void PID::BuildMappings(vector<UnpackerAcquConfig::hit_mapping_t>& hit_mappings,
                        vector<UnpackerAcquConfig::scaler_mapping_t>&) const
{
    for(const Element_t& element : elements)
    {
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC);

        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC);
    }
}



