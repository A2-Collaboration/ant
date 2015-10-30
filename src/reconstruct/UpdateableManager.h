#pragma once

#include "tree/TID.h"

#include <list>
#include <memory>

namespace ant {

class Updateable_traits;

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
    template<typename T>
    using shared_ptr_list = std::list< std::pair<std::size_t,std::shared_ptr<T>> >;

    std::list< std::pair< TID, shared_ptr_list<Updateable_traits> > > changePoints;
    const std::list< std::shared_ptr<Updateable_traits> > updateables;

    TID lastFlagsSeen;
};


}} // namespace ant::reconstruct
