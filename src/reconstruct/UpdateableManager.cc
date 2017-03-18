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


UpdateableManager::UpdateableManager(const std::list<std::shared_ptr<Updateable_traits>>& updateables_) :
    updateables(updateables_),
    lastFlagsSeen() // default ctor makes TID invalid
{

}

void UpdateableManager::UpdateParameters(const TID& currentPoint)
{
    // use last flags seen as some init flag
    if(lastFlagsSeen.IsInvalid()) {
        // ask each updateable for its items and build queue from it
        for(const shared_ptr<Updateable_traits>& updateable : updateables)
        {
            // tell starting TID
            updateable->UpdatedTIDFlags(currentPoint);

            // build queue from first call to Load
            for(auto item : updateable->GetLoaders()) {
                DoQueueLoad(currentPoint, item);
            }
        }
    }
    else if(currentPoint.Flags != lastFlagsSeen.Flags) {
        for(auto updateable : updateables)
            updateable->UpdatedTIDFlags(currentPoint);
    }
    lastFlagsSeen = currentPoint;

    // it might be that the current point lies far in the future
    // so calling Load more than once might be necessary
   while(!queue.empty() && queue.top().NextChangePoint <= currentPoint) {
       DoQueueLoad(queue.top().NextChangePoint, queue.top().Item);
       queue.pop();
   }
}

void UpdateableManager::DoQueueLoad(const TID& currPoint,
                                    Updateable_traits::Loader_t loader)
{
    TID nextChangePoint;
    loader(currPoint, nextChangePoint);
    if(nextChangePoint.IsInvalid())
        return;
    if(nextChangePoint <= currPoint) {
        LOG(WARNING) << "UpdateableItem returned NextChangePoint not pointing to the future";
        return;
    }
    queue.emplace(nextChangePoint, loader);
}
