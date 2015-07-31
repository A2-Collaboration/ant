#pragma once

#include <memory>

namespace ant {

class Event;
class TSlowControl;

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
    virtual bool HaveSlowControl() { return false; }
    virtual bool ReadNextEvent(Event& event, TSlowControl& slowControl) = 0;
};

}} //
