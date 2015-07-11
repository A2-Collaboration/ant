#include "UnpackerAcqu.h"
#include "detail/UnpackerAcqu_detail.h"

#include "tree/TDataRecord.h"
#include "base/Logger.h"

#include <stdexcept>

using namespace std;
using namespace ant;

UnpackerAcqu::UnpackerAcqu()
{

}


bool UnpackerAcqu::OpenFile(const std::string &filename)
{
  // this might also throw an exception if something
  // is strange with the file
  file = UnpackerAcquFileFormat::Get(filename, queue);
  // check if we were successful in finding a file
  if(file == nullptr)
    return false;

  LOG(INFO) << "Successfully opened " << filename;
  return true;
}

shared_ptr<TDataRecord> UnpackerAcqu::NextItem() noexcept
{
  // check if we need to replenish the queue
  if(queue.empty()) {
    file->FillEvents(queue);
    // still empty? Then the file is completely processed...
    if(queue.empty())
      return nullptr;
  }

  // std;:deque does not have a method to get and remove the element
  // we also convert the unique_ptr here to some shared_ptr
  const auto& element = shared_ptr<TDataRecord>(move(queue.front()));
  queue.pop_front();
  return element;
}


