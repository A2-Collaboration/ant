#include "SlowControlCreator.h"

using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;
using namespace std;

/**
 * @brief Get first value of vector. Returns NaN if empty.
 * @param cont Vector to access
 * @return Element 0, or NaN if vector is empty
 */
double SafeGet0(const std::vector<TKeyValue<double>>& cont) noexcept {
    return cont.empty() ? std::numeric_limits<double>::quiet_NaN() : cont[0].Value;
}

/**
 * @brief Unpack a TKetValue vector to linear vector of fixed size.
 *        Sets all elements to NaN and then fills elements according to the Key.
 *        Might throw std::out_of_range exception if Key >= size
 * @param source read from this
 * @param target write to this. size stays the same.
 */
void CopyVector(const std::vector<TKeyValue<double>>& source, std::vector<double>& target) {

    // set everything to nan
    target.assign(target.size(), std::numeric_limits<double>::quiet_NaN());

    // fill in values
    for(const auto& kv : source) {
            target.at(kv.Key) = kv.Value;
    }

}


void FillSlowControl(SlowControl& slc, const TSlowControl& value)
{

    switch (value.GetKey().Type) {

    case TSlowControl::Type_t::AcquScaler:

        if(value.Name == "Trigger/TotalLivetime") {
            slc.TotalLivetime() = SafeGet0(value.Payload_Float);
        }

        if(value.Name == "Trigger/FaradayCup") {
            slc.FaradayCup() = SafeGet0(value.Payload_Float);
        }

        if(value.Name == "Trigger/TaggerScalers") {
            try {
                CopyVector(value.Payload_Float, slc.TaggerScalers());
            } catch (const std::out_of_range&) {
                throw std::out_of_range("Tagger scaler channel number out of range");
            }
        }

        break;
    case TSlowControl::Type_t::EpicsOneShot:
        break;
    case TSlowControl::Type_t::EpicsScaler:
        break;
    case TSlowControl::Type_t::EpicsTimer:
        break;
    }
}


std::list<TSlowControl::Key> RequestedKeys(const SlowControl& slc)
{
    list<TSlowControl::Key> keys;

    if(slc.FaradayCup.isRequested())    { keys.emplace_back(TSlowControl::Key(TSlowControl::Type_t::AcquScaler, "Trigger/FaradayCup")); }
    if(slc.TotalLivetime.isRequested()) { keys.emplace_back(TSlowControl::Key(TSlowControl::Type_t::AcquScaler, "Trigger/TotalLivetime")); }
//    if(slc.TaggerScalers.isRequested()) { keys.emplace_back(TSlowControl::Key(TSlowControl::Type_t::AcquScaler), "Trigger/TaggerScalers"); }

    return keys;
}
