#include "UnpackerAcqu.h"
#include "UnpackerAcqu_detail.h"

#include "Logger.h"

using namespace std;
using namespace ant;

UnpackerAcqu::UnpackerAcqu()
{

}


bool UnpackerAcqu::OpenFile(const std::string &filename)
{
  // this might also throw an exception if something
  // is strange with the file
  file = UnpackerAcquFileFormat::Get(filename);
  // check if we were successful in finding a file
  if(file == nullptr)
    return false;



  LOG(INFO) << "Successfully opened " << filename;
  return true;
}

shared_ptr<TDataRecord> UnpackerAcqu::NextItem()
{
  return nullptr;
}

