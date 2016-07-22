#include "EPT.h"

#include "detail/EPT_2014_elements.h"

#include "tree/TID.h"

#include <limits>
#include <cassert>


using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;

const std::string EPT::ScalerName = "EPT_Scalers";

EPT::EPT(double beamEnergy, const std::vector<EPT::Element_t>& elements_init) :
    TaggerDetector_t(
        Detector_t::Type_t::EPT,
        beamEnergy
        ),
    elements(elements_init)
{
    assert(elements.size()==47);
}

void EPT::BuildMappings(
        vector<UnpackerAcquConfig::hit_mapping_t> &hit_mappings,
        vector<UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{

    vector<UnpackerAcquConfig::scaler_mapping_t::entry_t> scaler_entries;

    for(const Element_t& element : elements) {

        if(element.Ignored)
            continue;

        // TDC/scaler information is most important
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC
                                  );

        // build the scaler entries
        scaler_entries.emplace_back(element.Channel,
                                    element.Scaler);

        // ADC information are rarely present for the EPT
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC
                                  );
    }

    // map the scalers
    scaler_mappings.emplace_back(EPT::ScalerName, scaler_entries);
}





