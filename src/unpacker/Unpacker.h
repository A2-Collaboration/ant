#pragma once
/**
 * @page unpacker Unpacker
 * @tableofcontents
 *
 * The unpacker stage handles the conversion of raw data, for example AcquRoot
 * Mk2 data, to a stream of ant::T* classes derived from ant::TDataRecord, for
 * example ant::TDetectorRead.
 *
 * @section write_unpacker How to write a new unpacker?
 *
 * You need to implement the interfance ant::Unpacker::Module and then add the
 * instance to the list in the implementation of ant::Unpacker::Get.
 *
 * Please see ant::UnpackerA2Geant as a simple example based on ROOT tree
 * reading, and  ant::UnpackerAcqu for a more complicated example.
 *
 */

#include <string>
#include <memory>

namespace ant {

class TDataRecord;

/**
 * @brief The Unpacker class encapsulates the interface for unpacking raw data
 */
class Unpacker {

public:

    Unpacker() = delete;

    /**
     * @brief The Reader interface returns the unpacked data as a stream
     *
     * NextItem returns non-null pointers to instances of derived classes of TDataRecord,
     * usually TDetectorRead
     */
    class Reader {
    public:
        virtual ~Reader() = default;
        virtual std::unique_ptr<TDataRecord> NextItem() noexcept = 0;
        virtual double PercentDone() const = 0;
    };

    /**
     * @brief The Module interface is used by Unpacker::Get
     */
    class Module : public Reader {
        friend class Unpacker;
    protected:
        virtual bool OpenFile(const std::string& filename) = 0;
    };

    /**
     * @brief Get a unique unpacker instance for the given filename
     * @param filename the file to be examined for unpacking
     * @return pointer to the unpacker instance
     * @throw Exception if no or more than one unpacker for filename found
     */
    static std::unique_ptr<Reader> Get(const std::string &filename);

    /**
     * @brief The Exception class is thrown if an unexpected error during unpacking occurs
     */
    class Exception : public std::runtime_error {
        using std::runtime_error::runtime_error; // use base class constructor
    };
};

} // namespace ant
