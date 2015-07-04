#include "UnpackerAcqu.h"

#include "UnpackerAcqu_detail.h"
#include "RawFileReader.h"

#include "Logger.h"

#include "easylogging++.h"


#include "stl_helpers.h"
#include <iostream>

using namespace std;
using namespace ant;
using namespace ant::unpacker;

UnpackerAcqu::UnpackerAcqu()
{

}

//header_info_t GetHeaderInfo_Mk1(UInt_t* filebuffer) {
//  header_mk1_t* h = (header_mk1_t*)filebuffer;
//  header_info_t info;
//  info.Time = sanitize(h->fTime);
//  info.Description = sanitize(h->fDescription);
//  info.RunNote = sanitize(h->fRunNote);
//  info.OutFile = sanitize(h->fOutFile);
//  info.RunNumber = h->fRun;
//  info.IsMaybeMk2 =  *filebuffer == EHeadBuff;
//  info.IsSane  = true;
//  // check a few conditions for saneness
//  info.IsSane &= info.Time.length() == 24;
//  info.IsSane &= info.OutFile.length() != 0;
//  info.IsSane &= info.RunNumber >= 0;
//  return info;
//}

//header_info_t GetHeaderInfo_Mk2(vector<uint32_t> buffer) {
//  header_mk2_t* h = (header_mk2_t*)filebuffer;
//  header_info_t info;
//  info.Time = sanitize(h->fTime);
//  info.Description = sanitize(h->fDescription);
//  info.RunNote = sanitize(h->fRunNote);
//  info.OutFile = sanitize(h->fOutFile);
//  info.RunNumber = h->fRun;
//  info.IsMaybeMk2 =  *filebuffer == EHeadBuff;
//  info.IsSane  = true;
//  // check a few conditions for saneness
//  info.IsSane &= info.Time.length() == 24;
//  info.IsSane &= info.OutFile.length() != 0;
//  info.IsSane &= info.RunNumber >= 0;
//  return info;
//}


bool UnpackerAcqu::OpenFile(const std::string &filename)
{
  RawFileReader reader;
  reader.open(filename);
  vector<uint32_t> buffer(sizeof(acqu_t::AcquMk2Info_t)+1);
  reader.read(buffer.data(), buffer.size());
  if(buffer[0] != acqu_t::EHeadBuff)
    return false;

  // check for Mk2 revision

  acqu_t::AcquMk2Info_t* header =
      reinterpret_cast<acqu_t::AcquMk2Info_t*>(buffer.data()+1);

  //cout << std_ext::string_sanitize(header->fDescription);
  //cout << header->fMk2;


  LOG(INFO) << "Successfully opened " << filename;
  VLOG(2) << "Verbose!";
  return true;
}

shared_ptr<TDataRecord> UnpackerAcqu::NextItem()
{
  return nullptr;
}
