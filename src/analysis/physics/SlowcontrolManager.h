#pragma once

#include "tree/TSlowControl.h"
#include "analysis/data/Slowcontrol.h"

#include <map>
#include <queue>
#include <memory>
#include <list>

namespace ant {
namespace analysis {
namespace slowontrol {

class Manager {
protected:

    using buffer_t = std::queue<std::unique_ptr<const TSlowControl>>;
    std::map<TSlowControl::Key, buffer_t> slowcontrol;

    TID minimal_in_buffer;

    static TID min(const TID& a, const TID& b);

public:
    Manager() = default;

    void SetRequiredKeys(const std::list<TSlowControl::Key> keys);

    void ProcessSlowcontrol(std::unique_ptr<const TSlowControl> data);

    bool isComplete() const;

    TID FindMinimalTID() const;

    TID UpdateSlowcontrolData(data::Slowcontrol& slc);

};

}
}
}
