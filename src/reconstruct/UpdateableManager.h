#pragma once

#include <list>
#include <memory>

namespace ant {

class TID;
class Updateable_traits;

namespace reconstruct {

class UpdateableManager {
public:
    UpdateableManager(
            const TID& headerInfo_ID,
            const std::list< std::shared_ptr<Updateable_traits> >& updateables);

    void UpdateParameters(const TID& currentPoint);

private:
    template<typename T>
    using shared_ptr_list = std::list< std::shared_ptr<T> >;

    std::list< std::pair< TID, shared_ptr_list<Updateable_traits> > > changePoints;
};


}} // namespace ant::reconstruct
