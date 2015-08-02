#include "MemoryPool.h"

#include "tree/TDetectorRead.h"

using namespace std;
using namespace ant;

template<>
std::forward_list<std::unique_ptr<ant::TDetectorRead>> ant::MemoryPool<ant::TDetectorRead>::items{};
