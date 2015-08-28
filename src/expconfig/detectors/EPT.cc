#include "EPT.h"

#include "detail/EPT_2014_elements.h"

#include "tree/THeaderInfo.h"

#include <limits>
#include <cassert>


using namespace std;
using namespace ant;
using namespace ant::expconfig::detector;



bool EPT_2014::Matches(const THeaderInfo &headerInfo) const {
    return std_ext::time_between(headerInfo.Timestamp,
                                 "2014-07-01", "2014-12-31");
}

bool EPT::TryGetChannelFromPhoton(double photonEnergy, unsigned& channel) const {
    const double electronEnergy = BeamEnergy - photonEnergy;
    for(const Element_t& elem : elements) {
        const double lower = elem.ElectronEnergy - ElectronEnergyWidth/2;
        const double upper = elem.ElectronEnergy + ElectronEnergyWidth/2;
        if(electronEnergy >= lower && electronEnergy < upper) {
            channel = elem.Channel;
            return true;
        }
    }
    return false;
}

void EPT::BuildMappings(
        vector<UnpackerAcquConfig::hit_mapping_t> &hit_mappings,
        vector<UnpackerAcquConfig::scaler_mapping_t>& scaler_mappings) const
{
    for(const Element_t& element : elements) {
        // TDC/scaler information is most important
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Timing,
                                  element.Channel,
                                  element.TDC
                                  );
        scaler_mappings.emplace_back(Type,
                                     element.Channel,
                                     element.Scaler);

        // ADC information are rarely present for the EPT
        hit_mappings.emplace_back(Type,
                                  Channel_t::Type_t::Integral,
                                  element.Channel,
                                  element.ADC
                                  );


    }
}



