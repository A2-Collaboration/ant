#pragma once

#include "event_t.h"
#include "reader_flags_t.h"

namespace ant {
namespace analysis {
namespace input {

/**
 * @brief Abstract base class for data input modules
 *
 * Data input modules read MC/Detector data from somewhere.
 * Examples:
 *  * goat file reader
 *  * new ant data format reader
 *  * MCTrue Pluto reader
 */
class DataReader {
public:

    DataReader() = default;
    virtual ~DataReader() {}

    class Exception : public std::runtime_error {
      using std::runtime_error::runtime_error; // use base class constructor
    };

    virtual reader_flags_t GetFlags() const =0;
    virtual bool ReadNextEvent(event_t& event) =0;

    virtual double PercentDone() const =0;
};

}}} // namespace ant::analysis::input
