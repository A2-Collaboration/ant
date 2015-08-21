#pragma once

#include <memory>

namespace ant {
class TSlowControl;

namespace analysis {

namespace data {
class Event;
}

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

    DataReader() = default;
    virtual ~DataReader() {}

    class Exception : public std::runtime_error {
      using std::runtime_error::runtime_error; // use base class constructor
    };

    virtual bool IsSource() = 0;
    virtual bool ReadNextEvent(data::Event& event) = 0;
};

}}} // namespace ant::analysis::input
