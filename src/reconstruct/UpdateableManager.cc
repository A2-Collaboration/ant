#include "UpdateableManager.h"

#include "Reconstruct_traits.h"

#include "tree/TID.h"

#include "base/Logger.h"

#include <list>
#include <memory>
#include <map>
#include <stdexcept>

using namespace std;
using namespace ant;
using namespace ant::reconstruct;


UpdateableManager::UpdateableManager(const TID& startPoint,
        const std::list<std::shared_ptr<Updateable_traits> >& updateables_
        ) :
    updateables(updateables_),
    lastFlagsSeen(startPoint)
{

    // ask each updateable for its items and build queue from it
    for(const shared_ptr<Updateable_traits>& updateable : updateables)
    {
        // tell starting TID
        updateable->UpdatedTIDFlags(startPoint);

        // build queue from first call to Load
        for(auto item : updateable->GetItems()) {
            DoQueueLoad(startPoint, item);
        }
    }
}

void UpdateableManager::UpdateParameters(const TID& currentPoint)
{
    if(currentPoint.Flags != lastFlagsSeen.Flags) {
        for(auto updateable : updateables)
            updateable->UpdatedTIDFlags(currentPoint);
        lastFlagsSeen = currentPoint;
    }

    // it might be that the current point lies far in the future
    // so calling Load more than once might be necessary
   while(!queue.empty() && queue.top().ValidUntil < currentPoint) {
       DoQueueLoad(currentPoint, queue.top().Item);
       queue.pop();
   }
}

void UpdateableManager::DoQueueLoad(const TID& currPoint,
                                         Updateable_traits::UpdateableItemPtr item)
{
    TID validUntil;
    item->Load(currPoint, validUntil);
    if(validUntil.IsInvalid())
        return;
    if(validUntil < currPoint) {
        LOG(WARNING) << "UpdateableItem returned validity into the past";
        return;
    }
    queue.emplace(validUntil, item);
}
