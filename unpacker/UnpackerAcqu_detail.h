#ifndef UNPACKERACQU_DETAIL_H
#define UNPACKERACQU_DETAIL_H

#include <cstddef>
#include <cstdint>

// this header file contains copy-pasted code
// from the Acqu project below acqu_core/AcquRoot/src

// unused stuff was commented out but intentionally kept for reference
// all enum's were converted into constexpr inside the given namespace
// all types were replaced by definitions from <cstdint>

namespace ant {
namespace unpacker {
namespace acqu_t {

///////////////////
// EnumConst.h
///////////////////

//enum { EMk1 = 1, EMk2 };               // Acqu data formats
//enum { EMaxCmdList = 32 };             // max number different command lists
//enum { ELineSize = 1024 };             // length in characters of line of text
//enum { EKeyWordSize = 256 };           // length in characters of keyword
//enum { EErrNonFatal, EErrFatal };      // error level
//enum { ENOKeyWord, EKeyWord };         // Config file readin...line format
//enum { ENsem = 16 };                   // # of semaphores
//enum { EHostLen = 64 };	               // max length of host name (socket)
//enum { ESkDefPacket = 1024 };   // Default transfer length for stream socket
//enum { ESkLocal, ESkRemote };   // local/remote socket connect
//enum { ESkBacklog = 8 };        // max length pending connection queue (socket)
//enum { ESkInitBuff = 2 };       // # ints in initial handshake buffer (socket)
//enum { EMaxDataLength = 65536 };// max data buffer length
//enum { EMaxInputFiles = 1024 }; // size of input-file pointer buffer
//enum { EFalse, ETrue };         // Logic...should use kTRUE, kFALSE

// ACQU Mk1 data buffer header types and data delimiters
//enum{
constexpr uint32_t EHeadBuff = 0x10101010;      // header buffer (experimental parameters
constexpr uint32_t EDataBuff = 0x20202020;      // standard data buffer
constexpr uint32_t EEndBuff = 0x30303030;       // end-of-file buffer
constexpr uint32_t EKillBuff = 0x40404040;      // shut-down buffer
//  EPhysBuff = 0x50505050,      // reserved
//  EHeadPhysBuff = 0x60606060,  // reserved
constexpr uint32_t EMk2DataBuff = 0x70707070;   // Mk2 data buffer
constexpr uint32_t EEndEvent = 0xFFFFFFFF;      // end of event marker
constexpr uint32_t EBufferEnd = 0xFFFFFFFF;     // end of file marker
constexpr uint32_t EScalerBuffer = 0xFEFEFEFE;  // start of scaler read out
constexpr uint32_t EEPICSBuffer = 0xFDFDFDFD;   // start of EPICS read out
constexpr uint32_t EReadError = 0xEFEFEFEF;     // start of error block (hardware error)
//  ENullADC = -1,               // undefined ADC value
//  ENullHit = 0xFFFFFFFF,       // undefined hit index (end of hit buffer)
//  ENullStore = 0x8000,         // for multi-hit ADC handling
//  ENullFloat = -999999999      // optional null indicator
//};

// Constants for ROOT storage and analysis
//enum{ EMaxEventSize = 32768, EMaxName = 256 };

// Definitions for ADC setup
//enum{ EUndefinedADC = 0,          // ADC index not registered in analysis
//      EPatternADC = 0xffff,       // its a bit-pattern unit
//      EForeignADC = 0x3000,       // foreign data formats
//      EForeignScaler = 0x4000,    // ditto
//      EMultiADC = 0x10000,        // multi-hit ADC
//      EFlashADC = 0x20000         // flash ADC
//    };



////////////////////
// TA2Mk2Format.h
////////////////////

// Array sizes
constexpr size_t EMk2SizeTime = 32;
constexpr size_t EMk2SizeComment = 256;
constexpr size_t EMk2SizeFName = 128;
constexpr size_t EMk2SizeDesc = 256;

struct AcquMk2Info_t {	       	   // 1st part of header buffer
  uint32_t fMk2;                   // Mark as Mk2 data
  char fTime[EMk2SizeTime];	       // run start time (ascii)
  char fDescription[EMk2SizeDesc]; // description of experiment
  char fRunNote[EMk2SizeComment];  // particular run note
  char fOutFile[EMk2SizeFName];     // output file
  int32_t fRun;		              // run number
  int32_t fNModule;		      // total no. modules
  int32_t fNADCModule;	       	      // no. ADC modules
  int32_t fNScalerModule;	              // no. scaler modules
  int32_t fNADC;		              // no. ADC's read out
  int32_t fNScaler;	    	      // no. scalers readout
  int32_t fRecLen;		      // maximum buffer length = record len
};

struct ModuleInfoMk2_t {
  int32_t fModID;                 // Acqu ID # of module
  int32_t fModIndex;		// Module list index
  int32_t fModType;		// type of module, ADC, latch etc.
  int32_t fAmin;	       	        // 1st channel index
  int32_t fNChannel;	       	// # channels
  int32_t fNScChannel;            // # scaler channels (if ADC and Scaler)
  int32_t fBits;	       	        // significant bits from output word
};

struct ReadErrorMk2_t {
  uint32_t fHeader;	       	// error block header
  int32_t fModID;                 // hardware identifier
  int32_t fModIndex;	        // list index of module
  int32_t fErrCode;	       	// error code returned
  int32_t fTrailer;              // end of error block marker
};

// time was originally a time_t, now changed to UInt_t
// We found that time_t is 8-byte on a 64-bit machine and 4-byte on a
// 32-bit machine. Most front-end DAQ machines are 32-bit while
// most analysis/storage systems are 64-bit.
struct EpicsHeaderInfo_t {      //header for epics buffers and channels
  char name[32];              //Name of EPICs module
  uint32_t time;                  //Time of buffer
  int16_t index;                //index of module
  int16_t period; //0 = start only, +ve =period in counts, -ve =-1*period in ns
  int16_t id;                   //id of epics module (user set)
  int16_t nchan;                //no of EPICs channels in module
  int16_t len;                  //length of EPICS buffer
};

struct EpicsChannelInfo_t {      //header for channel info
  char pvname[32];            //Process variable name
  int16_t bytes;                //No of bytes for this channel
  int16_t nelem;                //No of elements in array
  int16_t type;                 //Type of element
};

//struct ScalerBlockMk2_t {
//  UInt_t fID;
//  UInt_t fCounts;
//}



}}} // namespace ant::unpacker::acqu_t

#endif // UNPACKERACQU_DETAIL_H
