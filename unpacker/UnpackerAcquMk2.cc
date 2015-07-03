#include "UnpackerAcquMk2.h"

#include "RawFileReader.h"


using namespace std;
using namespace ant;

UnpackerAcquMk2::UnpackerAcquMk2()
{

}

bool UnpackerAcquMk2::OpenFile(const std::string &filename)
{
  RawFileReader reader;
  reader.open(filename);
  vector<uint32_t> buffer(10);
  reader.read(buffer.data(), buffer.size());
  if(buffer[0] != 0x10101010)
    return false;



  return true;
}

shared_ptr<TDataRecord> UnpackerAcquMk2::NextItem()
{
  return nullptr;
}
