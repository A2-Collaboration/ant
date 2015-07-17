#include "catch.hpp"

#include "tmpfile_t.h"

#include "tree/TCalibrationData.h"
//#include "tree/TDataRecord.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

void dotest();

TEST_CASE("TCalibrationData: Write to TTree", "[tree]") {
  dotest();
}

void dotest() {
  tmpfile_t tmpfile;


  TFile f(tmpfile.filename.c_str(),"RECREATE");

  TTree* tree = new TTree("testtree","TCalibData Test Tree");
  ant::TCalibrationData* cdata = new ant::TCalibrationData(ant::TID(0,0,false),ant::TID(0,1,true));
//  ant::TCalibrationData* cdata = new ant::TCalibrationData();
}
