#include "ExpConfig.h"

#include "setups/Setup_2014_EtaPrime.h"

#include "unpacker/tree/THeaderInfo.h"
#include "unpacker/UnpackerAcqu.h"

#include "base/Logger.h"

#include <type_traits>
#include <list>

using namespace std;

namespace ant { // templates need explicit namespace

template<typename T>
unique_ptr<T> Get_(const THeaderInfo& header) {
  VLOG(9) << "Searching for config of type "
          << std_ext::getTypeAsString<T>();

  static_assert(is_base_of<ExpConfig::Base, T>::value, "T must be a base of ExpConfig::Base");

  // make a list of available configs
  std::list< std::unique_ptr<ExpConfig::Base> > modules;
  modules.push_back(std_ext::make_unique<expconfig::setup::Setup_2014_EtaPrime>());

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


unique_ptr<ExpConfig::Module> ExpConfig::Module::Get(const THeaderInfo& header)
{
  return Get_<ExpConfig::Module>(header);
}

unique_ptr<ExpConfig::Reconstruct> ExpConfig::Reconstruct::Get(const THeaderInfo& header)
{
  return Get_<ExpConfig::Reconstruct>(header);
}

template<>
unique_ptr< UnpackerAcquConfig >
ExpConfig::Unpacker<UnpackerAcquConfig>::Get(const THeaderInfo& header) {
  return Get_< UnpackerAcquConfig >(header);
}

} // namespace ant




