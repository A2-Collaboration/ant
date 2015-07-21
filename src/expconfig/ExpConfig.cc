#include "ExpConfig.h"

#include "setups/Setup_2014_EtaPrime.h"

#include "tree/THeaderInfo.h"
#include "unpacker/UnpackerAcqu.h"

#include "base/Logger.h"

#include <type_traits>
#include <list>

using namespace std;

namespace ant { // templates need explicit namespace

const std::list< std::shared_ptr<ExpConfig::Base> > ExpConfig::modules_available = {
  std::make_shared<expconfig::setup::Setup_2014_EtaPrime>()
};

template<typename T>
shared_ptr<T> ExpConfig::Get_(const THeaderInfo& header) {
  VLOG(9) << "Searching for config of type "
          << std_ext::getTypeAsString<T>();

  static_assert(is_base_of<Base, T>::value, "T must be a base of ExpConfig::Base");

  // make a copy of the list of available configs
  std::list< std::shared_ptr<Base> > modules = modules_available;

  // remove the config if the config says it does not match
  modules.remove_if([&header] (const shared_ptr<Base>& m) {
    return !m->Matches(header);
  });

  // check if something reasonable is left
  if(modules.empty()) {
    throw ExpConfig::Exception(std_ext::formatter()
                               << "No config found for header "
                               << header);
  }
  if(modules.size()>1) {
    throw ExpConfig::Exception(std_ext::formatter()
                               << "More than one config found for header "
                               << header);
  }

  // only one instance found, now try to cast it to the
  // requested type unique_ptr<T>
  // this only works if the found module is a derived class of the requested type

  const auto& ptr = dynamic_pointer_cast<T, Base>(modules.back());


  if(ptr==nullptr) {
    throw ExpConfig::Exception("Matched config does not fit to requested type");
  }

  // hand over the unique ptr
  return ptr;
}


shared_ptr<ExpConfig::Module> ExpConfig::Module::Get(const THeaderInfo& header)
{
  return Get_<ExpConfig::Module>(header);
}

shared_ptr<ExpConfig::Reconstruct> ExpConfig::Reconstruct::Get(const THeaderInfo& header)
{
  return Get_<ExpConfig::Reconstruct>(header);
}

template<>
shared_ptr< UnpackerAcquConfig >
ExpConfig::Unpacker<UnpackerAcquConfig>::Get(const THeaderInfo& header) {
  return Get_< UnpackerAcquConfig >(header);
}

} // namespace ant




