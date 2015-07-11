#include "ExpConfig.h"

#include "base/Logger.h"

#include "unpacker/tree/THeaderInfo.h"
#include "unpacker/UnpackerAcqu.h"

#include "setups/Setup_2014_EtaPrime.h"

#include <type_traits>
#include <list>
#include <cxxabi.h>
#include <stdexcept>

using namespace std;
using namespace ant;

bool Channel_t::IsIntegral(const Channel_t::Type_t& t) {
  switch(t) {
  case Type_t::Integral:
  case Type_t::IntegralShort:
  case Type_t::IntegralAlternate:
  case Type_t::IntegralShortAlternate:
    return true;
  default:
    return false;
  }
}

string Channel_t::ToString(const Type_t& type)
{
  switch(type) {
  case Channel_t::Type_t::BitPattern:
    return "BitPattern";
  case Channel_t::Type_t::Counter:
    return "Counter";
  case Channel_t::Type_t::Integral:
    return "Integral";
  case Channel_t::Type_t::IntegralAlternate:
    return "IntegralAlternate";
  case Channel_t::Type_t::IntegralShort:
    return "IntegralShort";
  case Channel_t::Type_t::IntegralShortAlternate:
    return "IntegralShortAlternate";
  case Channel_t::Type_t::Scaler:
    return "Scaler";
  case Channel_t::Type_t::Timing:
    return "Timing";
  }
  throw runtime_error("Not implemented");
}

string Detector_t::ToString(const Type_t &type)
{
  switch(type) {
  case Detector_t::Type_t::CB :
    return "CB";
  case Detector_t::Type_t::Cherenkov:
    return "Cherenkov";
  case Detector_t::Type_t::EPT:
    return "EPT";
  case Detector_t::Type_t::Moeller:
    return "Moeller";
  case Detector_t::Type_t::MWPC0:
    return "MWPC0";
  case Detector_t::Type_t::MWPC1:
    return "MWPC1";
  case Detector_t::Type_t::PID:
    return "PID";
  case Detector_t::Type_t::Tagger:
    return "Tagger";
  case Detector_t::Type_t::TAPS:
    return "TAPS";
  case Detector_t::Type_t::TAPSVeto:
    return "TAPSVeto";
  case Detector_t::Type_t::Trigger:
    return "Trigger";
  }
  throw runtime_error("Not implemented");
}


namespace ant { // templates need explicit namespace

template<typename T>
unique_ptr<T> Get_(const THeaderInfo& header) {
  VLOG(9) << "Searching for config of type "
          << abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, nullptr);

  static_assert(is_base_of<ExpConfig::Base, T>::value, "T must be a base of ExpConfig::Base");

  // make a list of available configs
  std::list< std::unique_ptr<ExpConfig::Base> > modules;

  modules.emplace_back(new expconfig::setup::Setup_2014_EtaPrime());

  // remove the config if the config says it does not match
  modules.remove_if([&header] (const unique_ptr<ExpConfig::Base>& m) {
    return !m->Matches(header);
  });

  // check if something reasonable is left
  if(modules.empty()) {
    throw ExpConfig::Exception("No config found for header "+header.Description);
  }
  if(modules.size()>1) {
    throw ExpConfig::Exception("More than one config found for header "+header.Description);
  }

  // only one instance found, now try to cast it to the
  // requested type unique_ptr<T>
  // this only works if the found module is a derived class of the requested type

  T* ptr = dynamic_cast<T*>(modules.back().release());

  if(ptr==nullptr) {
    throw ExpConfig::Exception("Matched config does not fit to requested type");
  }

  // hand over the unique ptr
  return unique_ptr<T>(ptr);
}


unique_ptr<ExpConfig::Module> ExpConfig::Get(const THeaderInfo& header)
{
  return Get_<ExpConfig::Module>(header);
}

template<>
unique_ptr< UnpackerAcquConfig >
ExpConfig::Unpacker<UnpackerAcquConfig>::Get(const THeaderInfo& header) {
  return Get_< UnpackerAcquConfig >(header);
}





} // namespace ant




