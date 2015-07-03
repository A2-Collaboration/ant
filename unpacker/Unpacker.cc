#include "Unpacker.h"
#include "UnpackerAcquMk2.h"

#include <algorithm>

#include <iostream>

using namespace std;
using namespace ant;

Unpacker::Unpacker()
{
  modules.emplace_back(new UnpackerAcquMk2());
}

unique_ptr<Unpacker::Module> Unpacker::Get(const string& filename)
{
  Unpacker unpacker;
  unpacker.modules.remove_if([&filename] (const unique_ptr<Unpacker::Module>& m) {
    return !m->OpenFile(filename);
  });
  if(unpacker.modules.empty()) {
    throw Exception("No suitable unpacker found for file "+filename);
  }
  if(unpacker.modules.size()>1) {
    throw Exception("More than one unpacker found for file "+filename);
  }
  // hand over the unique ptr
  return std::move(unpacker.modules.back());
}
