#include "ExpConfig.h"

#include "Logger.h"

#include "THeaderInfo.h"

#include "UnpackerAcqu.h"

using namespace std;
using namespace ant;



#include <list>


shared_ptr<ExpConfig::Module> ExpConfig::Get(const THeaderInfo& header)
{
//  // make a list of available configs
//  std::list< std::unique_ptr<Module> > modules;
//  //modules.emplace_back(new UnpackerAcqu());

//  // remove the config if the config says it does not match
//  modules.remove_if([&header] (const unique_ptr<Module>& m) {
//    return !m->Matches(header);
//  });

//  // check if something reasonable is left
//  if(modules.empty()) {
//    throw Exception("No config found for header "+header.Description);
//  }
//  if(modules.size()>1) {
//    throw Exception("More than one config found for header "+header.Description);
//  }

//  // hand over the unique ptr
//  return std::move(modules.back());
  return nullptr;
}

namespace ant { // templates need explicit namespace

namespace config {

namespace detector {

class TAPS {

};

}

namespace beamtime {

class EtaPrime :
    public ExpConfig::Module,
    public ExpConfig::Unpacker<UnpackerAcqu>
{


};

} // namespace beamtime
} // namespace configs


template<>
shared_ptr< ExpConfig::Unpacker<UnpackerAcqu> >
ExpConfig::Unpacker<UnpackerAcqu>::Get(const THeaderInfo& headerInfo) {

  VLOG(9) << "Hello";

  return nullptr;
}

} // namespace ant
