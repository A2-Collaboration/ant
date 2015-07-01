#ifndef DATAREADER_H
#define DATAREADER_H

#include "Event.h"
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
};

}

}
#endif
