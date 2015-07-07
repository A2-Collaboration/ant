#ifndef DATAREADER_H
#define DATAREADER_H

#include "data/Event.h"
#include <memory>

namespace ant {
namespace input {

/**
 * @brief Abstract base class for data input modules
 *
 * Data input modules read MC/Detector data from somewhere.
 * Examples:
 *  * goat file reader
 *  * new ant data foramt reader
 */
class DataReader {
public:
    virtual ~DataReader() = default;
    virtual std::shared_ptr<Event> ReadNextEvent() =0;
    virtual bool hasData() const =0;


    virtual long long EventsRead() const =0;

    /**
     * @brief Get total number of events available
     * @return number, -1 if unknown/does not apply
     */
    virtual long long TotalEvents() const;
};

}

}
#endif
