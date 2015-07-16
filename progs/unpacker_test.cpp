
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "tree/THeaderInfo.h"
#include "tree/TUnpackerMessage.h"
#include "tree/TSlowControl.h"
#include "tree/TDetectorRead.h"

#include "unpacker/RawFileReader.h"

#include "reconstruct/Reconstruct.h"

#include "base/Logger.h"
#include "base/Format.h"

#include <chrono>

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace ant;

int main(int argc, char* argv[]) {
  SetupLogger(argc, argv);

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();

  //auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_9227.dat");
  //auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_7892.dat");
  auto unpacker = Unpacker::Get("scratch/CBTaggTAPS_5711.dat.xz");
  //auto unpacker = Unpacker::Get("scratch/oneevent-small.dat");
  unsigned nReads = 0;

  // write some stuff to a ROOT tree
  TFile* file = new TFile("unpacker_test_out.root","RECREATE");

  TTree* treeHeaderInfo = new TTree("treeHeaderInfo", "treeHeaderInfo");
  THeaderInfo* HeaderInfo = new THeaderInfo();
  treeHeaderInfo->Branch("THeaderInfo", "ant::THeaderInfo", &HeaderInfo);

  TTree* treeUnpackerMessage = new TTree("treeUnpackerMessage", "treeUnpackerMessage");
  TUnpackerMessage* UnpackerMessage = new TUnpackerMessage();
  treeUnpackerMessage->Branch("UnpackerMessage", "ant::TUnpackerMessage", &UnpackerMessage);

  TTree* treeSlowControl = new TTree("treeSlowControl", "treeSlowControl");
  TSlowControl* SlowControl = new TSlowControl();
  treeSlowControl->Branch("SlowControl", "ant::TSlowControl", &SlowControl);

  TTree* treeDetectorRead = new TTree("treeDetectorRead", "treeDetectorRead");
  TDetectorRead* DetectorRead = new TDetectorRead();
  treeDetectorRead->Branch("DetectorRead", "ant::TDetectorRead", &DetectorRead);


  unique_ptr<Reconstruct> reconstruct;

  while(auto item = unpacker->NextItem()) {
    if(nReads>30000)
      break;

    //cout << *item << endl;
    HeaderInfo = dynamic_cast<THeaderInfo*>(item.get());
    if(HeaderInfo != nullptr) {
      treeHeaderInfo->Fill();
      reconstruct = std_ext::make_unique<Reconstruct>(*HeaderInfo);
      continue;
    }
    UnpackerMessage = dynamic_cast<TUnpackerMessage*>(item.get());
    if(UnpackerMessage != nullptr) {
      treeUnpackerMessage->Fill();
      continue;
    }
    SlowControl = dynamic_cast<TSlowControl*>(item.get());
    if(SlowControl != nullptr) {
      treeSlowControl->Fill();
      continue;
    }
    DetectorRead = dynamic_cast<TDetectorRead*>(item.get());
    if(DetectorRead != nullptr) {

      if(reconstruct)
        reconstruct->DoReconstruct(*DetectorRead);

      treeDetectorRead->Fill();

      nReads++;
      if(nReads % 10000 == 0) {
        VLOG(5) << "Unpacked " << nReads << " detector reads";
      }
      continue;
    }
  }

  file->Write();
  file->Close();
  delete file;

  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  cout << "Analyzed " << nReads << " reads, speed " << nReads/elapsed_seconds.count() << " Reads/s" << endl;

  return EXIT_SUCCESS;
}
