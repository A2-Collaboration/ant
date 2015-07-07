#ifndef UNPACKERACQU_LEGACY_H
#define UNPACKERACQU_LEGACY_H

#include <cstdint>
#include <map>
#include <string>


// the following is copy-pasted code
// from the Acqu project below acqu_core/AcquRoot/src

// unused stuff was commented out but intentionally kept for reference
// all enum's were converted into constexpr inside the given namespace
// all types were replaced by definitions from <cstdint>

namespace ant {
namespace unpacker {
namespace acqu {

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
constexpr int32_t EUndefinedADC = 0;          // ADC index not registered in analysis
constexpr int32_t EPatternADC = 0xffff;       // its a bit-pattern unit
constexpr int32_t EForeignADC = 0x3000;       // foreign data formats
constexpr int32_t EForeignScaler = 0x4000;    // ditto
constexpr int32_t EMultiADC = 0x10000;        // multi-hit ADC
constexpr int32_t EFlashADC = 0x20000;        // flash ADC

/////////////////////
// DataFormats.h
/////////////////////

//	ACQU experimental header struct
//	formerly called m_dev_header

struct AcquExptHeader_t {
  char time[26];		// time in ascii
  char exp_desc[133];		// experiment header
  char run_note[133];		// particular run note
  char out_file[40];		// output file
  uint16_t run;			// run number
  uint16_t n_slave;		// no. of slave VME crates
  uint16_t n_module;		// total no. modules
  uint16_t n_vme;	      	// act. no. VME modules
  uint16_t n_camac;		// act. no. CAMAC modules
  uint16_t n_fastb;		// act. no. FASTBUS modules
  uint16_t n_spect;		// act. no readout subaddrs
  uint16_t n_scaler;		// no. of scalers readout
  uint16_t C_IRQ;		// no. CAMAC IRQ readout addresses
  uint16_t C_IRQ_S;		// no. CAMAC scaler readout addresses
  uint16_t F_IRQ;		// no. FASTBUS modules for IRQ readout
  uint16_t F_IRQ_S;		// no. FASTBUS modules for scaler readout
  uint16_t len_buff;		// length of data buff(byte)
};

struct AcquExptInfo_t {	       	// 1st part of header buffer
  char fTime[26];		// run start time (ascii)
  char fDescription[133];	// description of experiment
  char fRunNote[133];		// particular run note
  char fOutFile[40];		// output file
  uint16_t fRun;		// run number
  uint16_t fNslave;		// no. of slave VME crates
  uint16_t fNmodule;		// total no. modules
  uint16_t fNvme;	       	// no. VME modules
  uint16_t fNcamac;		// no. CAMAC modules
  uint16_t fNfastb;		// no. FASTBUS modules
  uint16_t fNspect;		// no. ADC's read out
  uint16_t fNscaler;		// no. scalers readout
  uint16_t fCIrq;	       	// no. CAMAC ADC module readouts
  uint16_t fCIrqS;		// no. CAMAC scaler module readouts
  uint16_t fFIrq;	       	// no. FASTBUS ADC module readouts
  uint16_t fFIrqS;		// no. FASTBUS scaler module readouts
  uint16_t fRecLen;		// maximum buffer length = record len
};

//	This struct specifies the module from which
//	an ADC or scaler has been read
//	formerly called P_spect

struct ADCInfo_t {
  uint16_t fModIndex;		// index in ModHeader array
  uint16_t fModSubAddr;		// subaddress in module
};
struct Pspect_t {
  uint16_t sp_index;		// index of module ID structure
  uint16_t a;			// module subaddress
};

//	Header struct of a hardware module...
//	Information about a particular piece of hardware
//	formerly called mod_header

struct ModHeader_t{
  char name[20];		// name of module
  uint16_t vic_crate;		// VME crate no. on VIC bus
  uint16_t bus_type;		// index vme, camac, fastbus etc.
  uint16_t mod_type;		// type of module, ADC, latch etc.
  uint16_t branch;		// branch address
  uint16_t crate;	       	// crate address
  uint16_t station;		// crate station address
  uint16_t a_min;	       	// min subaddress
  uint16_t a_max;	       	// max subaddress
  uint16_t bits;	       	// max no. bits from output
};

struct ModuleInfo_t {
  char fName[20];		// name of module
  uint16_t fVicCrate;		// VME crate no. on VIC bus
  uint16_t fBusType;		// index vme, camac, fastbus etc
  uint16_t fModType;		// type of module, ADC, latch etc.
  uint16_t fBranch;		// branch address (geographical)
  uint16_t fCrate;		// crate address (geographical)
  uint16_t fStation;		// crate station address (geographical)
  uint16_t fAmin;	       	// 1st subaddress
  uint16_t fAmax;	       	// last subaddress
  uint16_t fBits;	       	// max no. bits from output word
};

//	Structure of an error block...written if the front-end readout
//	detects an error
//	formerly called read_error

struct ReadError_t {
  uint32_t fHeader;	       	// error block header
  uint16_t fBus;	       	// bus type
  uint16_t fCrate;		// crate no. (geographical)
  uint16_t fStation;		// station no. (geographical)
  uint16_t fSubAddr;		// subaddress (geographical)
  uint16_t fCode;	       	// error code returned
  uint16_t fAlign;		// to keep long word boundary alignment
};

//	For converting foreign data formats to ACQU

struct AcquBlock_t{
  uint16_t id;			// adc index
  uint16_t adc;			// adc value
};
//struct InvAcquBlock_t{
//  uint16_t adc;			// adc value
//  uint16_t id;			// adc index
//};

//---------------------------------------------------------------------------
//	Lund Data Format
//	H.Ruijter data acquisition system
//--------------------------------------------------------------------------

//#define		LUND_REC_LEN	16384		// default record length
//#define		ADC_BLOCK	0xFDFDFDFD	// ADC readout event header
//#define		SCALER_BLOCK	0xFEFEFEFE	// Scaler readout event header
//#define		MAX_LUND_ADC	1024		// Maximum no. ADC's read

//	Struct of header at start of each event in data buffer

//struct LundEventHeader_t{
//  uint32_t header;	        // FDFDFDFD for adc; FEFEFEFE for scaler
//  uint32_t res0;		        // "reserved"
//  uint32_t n;		        // no. adc's or scalers
//  uint32_t res1;		        // "reserved"
//  uint32_t event;		        // event number
//};

//	Struct of start of header buffer
//	Differs slightly from standard ACQU

//struct LundExptHeader_t{
//  char time[28];	       	// time in ascii
//  char exp_desc[128];	       	// experiment description
//  char run_note[128];	       	// run description
//  char out_file[40];	       	// data file name
//  uint32_t run;		       	// run number
//};

//	End of Lund-specific stuff

//---------------------------------------------------------------------------
//	TAPS data format specification
//---------------------------------------------------------------------------

//	Comes at the start of every 8k buffer
//	Note there is a difference between Experimental data
//	and Geant-simulation data

//struct TapsExptHeader_t{
//  uint16_t lbuffer;	       	// buffer length..exp data
//  uint16_t nbuff;      		// buffer counter..exp data
//  uint16_t id;	       		// ID = 4
//  uint16_t polarisation; 	// beam polarisation info
//  uint16_t lgbuffer;	       	// buffer length..geant data
//  uint32_t nbuffer;	      	// buffer number
//  uint16_t reserved[16];       	// 32 bytes presently unused
//};

//	At the start of every event in the buffer

//struct TapsEventHeader_t{
//  int16_t levent;	       	// event length
//  uint16_t fragment;	       	// event fragmented?
//  uint16_t id;		       	// ID = 4
//  uint16_t info;				// scaler and polarisation info
//};

//	Sub-event specification..at the start of every sub-event within
//	an event....differences between experimental and simulation data

//struct TapsSubEvent_t{	       	// experimental header
//  uint16_t levent;	       	// length in bytes incl header
//  uint16_t id;		       	// subevent type id
//};

//struct GTapsSubEvent_t{	        // simulation header
//  uint16_t id;		       	// subevent type id
//  uint16_t levent;	       	// length shorts after header
//};

//	Standard TAPS grouping of index and 3 ADC's

//struct TapsBlock_t{	       	// 1 of these for each hit det.
//  uint16_t id;		       	// id number
//  uint16_t adc[3];	       	// TDC, QDC, QDC
//};

//	FERA readout TAPS grouping of index and 2 ADC's

//struct TapsFera_t{	       	// 1 of these for each hit det.
//  uint16_t id;		       	// id number
//  uint16_t adc[2];	       	// TDC, QDC
//};

//	Veto detector grouping of index and ADC

//struct TapsVeto_t{	       	// 1 for every hit veto counter
//  uint16_t id;		       	// id number
//  uint16_t adc;		       	// charge output
//};

//	Sub-event id = 25. Mainz tagger stuff

//struct TapsTagger_t{
//  uint16_t nevent;	       	// tagger event number
//  uint16_t levent;	       	// tagger event length
//  uint16_t sevent;	       	// sub-sub event spec.
//};
//struct TapsTagger1_t{	       	// sevent = 1
//  uint16_t id;		       	// channel id
//  uint16_t time;	       	// tdc output
//};
//struct TapsTagger2_t{	       	// sevent = 2
//  uint16_t taggeff;	       	// taggeff photon detector QDC
//};
//struct TapsTagger3_t{	       	// sevent = 3
//  uint32_t scalers[352];	       	// tagger scalers
//};
//struct TapsTagger4_t{	       	// sevent = 4
//  uint32_t scalers[12];	       	// tagger camac scalers
//};
//struct TapsTagger5_t{	       	// sevent = 5
//  uint16_t latch[24];	       	// tagger latches
//};


//---------------------------------------------------------------------------
//	Data format specification for Giessen test data acquisition
//	system. NB data is a 2-byte unsigned int...here referred to
//	as a word.
//	The header struct comes at the start of each buffer. Individual
//	events are fixed length and have no header
//---------------------------------------------------------------------------

//struct GtestEventHeader_t{
//  uint16_t id;		       	// buffer identification
//  uint16_t nbuffer;	       	// buffer number
//  uint16_t lbuffer;	       	// used no. of words in buffer
//  uint16_t levent;	       	// no. words in an event
//  uint32_t time;		       	// time (seconds since 1/1/70)
//  uint32_t reserved;	       	// not used
//};


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

//////////////////////////////
// Adapted from ModuleIndex.h
//////////////////////////////

constexpr uint32_t EDAQ_Undef=0;
constexpr uint32_t EDAQ_ADC=0x1;
constexpr uint32_t EDAQ_Scaler=0x2;
constexpr uint32_t EDAQ_SlowCtrl=0x4;
constexpr uint32_t EDAQ_Ctrl=0x8;
constexpr uint32_t EDAQ_SCtrl=0x10;
constexpr uint32_t EDAQ_IRQ=0x20;
constexpr uint32_t EDAQ_VADC=0x40;
constexpr uint32_t EDAQ_VScaler=0x80;
//constexpr uint32_t EDAQ_ADC_Scaler=0x3;
constexpr uint32_t EDAQ_Epics=0x100;

// add IDs as you encounter them in files
// do not cargo-cult totally unused IDs

static const std::map<uint32_t, std::string> ModuleIDToString = {
  // Controllers
//  {"CAEN_V2718", ECAEN_V2718},
  {0xf001, "KPH_I686"},
  // VMEbus
//  {"VMEbus",     EVMEbus},
  {0x0001, "CBD_8210"},
//  {"HSM_8170",   EHSM_8170},
//  {"VIC_8250",   EVIC_8250},
//  {"SLV_8250",   ESLV_8250},
  {0x0100, "LRS_1131"},
  {0x0400, "KPH_SMI"},
  {0x0401, "KPH_VITEC"},
//  {"CAEN_1465",  ECAEN_1465},
//  {"CAEN_1488",  ECAEN_1488},
  {0x0302, "CAEN_V792"},
//  {"CAEN_V775",  ECAEN_V775},
  {0x304, "CAEN_V874"},
  {0x305, "CAEN_V1190"},
  {0x306, "CAEN_V965"},
//  {"CAEN_SY1527",ECAEN_SY1527},
//  {"CAEN_208",   ECAEN_208},
//  {"INFN_56001", EINFN_56001},
//  {"FIADC_64",   EFIADC_64},
  {0x601, "CATCH_TDC"},
//  {"CATCH_SCA",  ECATCH_SCA},
  {0x0603, "GeSiCA"},
  {0x604, "CATCH_TCS"},
  {0x0700, "GSI_VUPROM"},
  {0x0701, "GSI_VUPROMT"},
  {0x0702, "GSI_VUPROM_Scaler"},
  {0x0703, "GSI_VUPROM_Moeller"},
  {0x0704, "GSI_VUPROM_Pattern"},
//  {"SIS_3820",   ESIS_3820},
//  // CAMAC
  {0x1fff, "CAMAC"},
//  {"VUCAM",      EKPH_VUCAM},
//  {"LRS_4508",   ELRS_4508},
//  {"LRS_4516",   ELRS_4516},
//  {"LRS_4434",   ELRS_4434},
//  {"LRS_2341",   ELRS_2341},
//  {"LRS_2373",   ELRS_2373},
//  {"LRS_3420",   ELRS_3420},
//  {"LRS_2323",   ELRS_2323},
//  {"LRS_2891",   ELRS_2891},
//  {"LRS_2132",   ELRS_2132},
//  {"LRS_4413",   ELRS_4413},
//  {"LRS_4413SPILL",   ELRS_4413SPILL},
//  {"LRS_4413DAQEnable",   ELRS_4413DAQEnable},
//  {"GSI_4800",   EGSI_4800},


//  // FASTBUS
//  {"FASTBUS",    EFASTBUS},
//  {"LRS_1821",   ELRS_1821},
//  {"LRS_1875",   ELRS_1875},
//  {"LRS_1872",   ELRS_1872},
//  {"LRS_1885",   ELRS_1800},
//  {"LRS_1877",   ELRS_1877},
//  {"LRS_1879",   ELRS_1879},
//  {"LRS_1881",   ELRS_1881},
//  {"PH_10C6",    EPH_10C6},
//  {"PH_10C2",    EPH_10C2},
//  {"STR_200",    ESTR_200},
//  {"STR_136",    ESTR_136},
//  // I2C
//  // PCI
//  {"Virtual",    EDAQ_Virtual},
//  {"EPICS",      EDAQ_Epics},
};

}}} // namespace ant::unpacker::acqu

#endif // UNPACKERACQU_LEGACY_H
