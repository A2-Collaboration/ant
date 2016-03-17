#include "UnpackerAcqu.h"
#include "detail/UnpackerAcqu_detail.h"

#include "tree/TEvent.h"
#include "base/Logger.h"

#include <stdexcept>

using namespace std;
using namespace ant;

UnpackerAcqu::UnpackerAcqu() {}
UnpackerAcqu::~UnpackerAcqu() {}

double UnpackerAcqu::PercentDone() const
{
    return file->PercentDone();
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

TEvent UnpackerAcqu::NextEvent() noexcept
{
    // check if we need to replenish the queue
    if(queue.empty()) {
        file->FillEvents(queue);
        // still empty? Then the file is completely processed...
        if(queue.empty())
            return {};
    }

    // std;:deque does not have a method to get and remove the element
    auto element = move(queue.front());
    queue.pop_front();
    return element;
}





