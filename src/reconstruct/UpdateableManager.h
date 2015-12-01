#pragma once

#include "tree/TID.h"

#include "Reconstruct_traits.h"

#include <queue>
#include <list>
#include <memory>

namespace ant {


namespace reconstruct {

class UpdateableManager {
public:
    /**
     * @brief UpdateableManager initializes the manager
     * @param startPoint is the minimum time point from which parameters are needed
     * @param updateables list of updateable items to be managed
     */
    UpdateableManager(
            const TID& startPoint,
            const std::list< std::shared_ptr<Updateable_traits> >& updateables_);

    /**
     * @brief UpdateParameters make the managed items ready for given currentPoint
     * @param currentPoint the time point
     */
    void UpdateParameters(const TID& currentPoint);

private:
    struct queue_item_t {
        TID ValidUntil;
        Updateable_traits::UpdateableItemPtr Item;
        queue_item_t(const TID& validUntil,
                     const Updateable_traits::UpdateableItemPtr& item) :
            ValidUntil(validUntil),
            Item(item)
        {}
        bool operator<(const queue_item_t& other) const {
            // invert ordering such that item with earliest TID
            // comes first in priority queue
            return other.ValidUntil < ValidUntil;
        }
    };

    std::priority_queue<queue_item_t> queue;

    std::list< std::shared_ptr<Updateable_traits> > updateables;
    TID lastFlagsSeen;

    void DoQueueLoad(const TID& currPoint,
                          Updateable_traits::UpdateableItemPtr item);
};


}} // namespace ant::reconstruct
