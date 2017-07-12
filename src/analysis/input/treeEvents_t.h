#pragma once

#include "tree/TEvent.h"
#include "base/WrapTTree.h"

namespace ant {
namespace analysis {
namespace input {

struct treeEvents_t : WrapTTree {
    ADD_BRANCH_T(TEvent, data)
};

}}}