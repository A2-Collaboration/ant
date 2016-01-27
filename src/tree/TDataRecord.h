#pragma once

/** @page tree Tree
 *
 * This part of the project contains the classes used to read/write data using
 * ROOT trees.
 *
 * @section unpackerstream Unpacker stream classes
 *
 * The classes
 *   - ant::TDetectorRead,
 *   - ant::TUnpackerMessage,
 *   - ant::THeaderInfo,
 *   - ant::TEvent and
 *   - ant::TSlowControl
 *
 * all derive from ant::TDataRecord and are used by the unpackers to return the
 * unpacked raw data as some well formatted data stream. The ant::TDetectorRead
 * items are then processed by the reconstruct stage.
 *
 * All those data stream classes are stored in separate trees by the
 * ant::tree::UnpackerWriter and can be read in again by the
 * ant::tree::UnpackerReader.
 *
 * @section calibration Calibration
 *
 * The class ant::TCalibrationData represents one saved calibration iteration.
 *
 */

#include "TID.h"

#include "Rtypes.h"

#include <iomanip>
#include <ctime>

namespace ant {

#ifndef __CINT__
struct TDataRecord : printable_traits
#else
struct TDataRecord
#endif
{
    TDataRecord() : ID() {}
    TDataRecord(const TID& id) : ID(id) {}
    virtual ~TDataRecord() {}

    TID ID;

#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override {
        return s << "TDataRecord ID=" << ID;
    }
#endif

    ClassDef(TDataRecord, ANT_UNPACKER_ROOT_VERSION)

}; // TDataRecord


} // namespace ant
