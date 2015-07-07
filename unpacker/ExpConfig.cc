#include "ExpConfig.h"

#include "Logger.h"

#include "THeaderInfo.h"

//#include "UnpackerAcqu.h"

using namespace std;
using namespace ant;


#include <type_traits>
#include <list>


//template<>
//unique_ptr<ExpConfig::Module> Get(const THeaderInfo& headerInfo);



namespace ant { // templates need explicit namespace

class UnpackerAcqu;

namespace config {

namespace detector {

struct Detector {
  const Detector_t Type;
protected:
  Detector(const Detector_t& type) :
    Type(type) {}
  Detector(const Detector&) = delete; // forbid copy
  virtual ~Detector() = default;
};

struct CB : public Detector {
  CB() : Detector(Detector_t::CB) {}
};

struct TAPS : public Detector {
  TAPS() : Detector(Detector_t::TAPS) {}
};

} // namespace detector

namespace setup {

} // namespace setup

namespace beamtime {

class EtaPrime :
    public ExpConfig::Module,
    public ExpConfig::Unpacker<UnpackerAcqu>
{

  // Base interface
public:
  bool Matches(const THeaderInfo &header) const {
    return true;
  }
};


} // namespace beamtime

} // namespace config




template<typename T>
unique_ptr<T> Get_(const THeaderInfo& header) {
  static_assert(is_base_of<ExpConfig::Base, T>::value, "T must be a base of ExpConfig::Base");

  // make a list of available configs
  std::list< std::unique_ptr<ExpConfig::Base> > modules;

  modules.emplace_back(new config::beamtime::EtaPrime());

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
    throw ExpConfig::Exception("Matched config does not fit to type");
  }

  // hand over the unique ptr
  return unique_ptr<T>(ptr);
}


unique_ptr<ExpConfig::Module> ExpConfig::Get(const THeaderInfo& header)
{
  return Get_<ExpConfig::Module>(header);
}

template<>
unique_ptr< ExpConfig::Unpacker<UnpackerAcqu> >
ExpConfig::Unpacker<UnpackerAcqu>::Get(const THeaderInfo& header) {

  VLOG(9) << "Hello";

  return Get_< ExpConfig::Unpacker<UnpackerAcqu> >(header);
}


} // namespace ant




