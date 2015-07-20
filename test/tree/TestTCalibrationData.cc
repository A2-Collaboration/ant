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
  tree->Branch("cdata",cdata);


  cdata->Data.push_back(ant::TCalibrationEntry(1,2.1));
  tree->Fill();
  cdata->Data.push_back(ant::TCalibrationEntry(2,3.2));
  tree->Fill();


  cout << cdata << endl;

  f.Write();
  f.Close();

  tree = nullptr;

  TFile f2(tmpfile.filename.c_str(),"READ");
  REQUIRE(f2.IsOpen());

  f2.GetObject("testtree",tree);

  REQUIRE(tree!=nullptr);

  ant::TCalibrationData* readcdata = nullptr;
  tree->SetBranchAddress("cdata",&readcdata);

  REQUIRE(tree->GetEntries() == 2);

  tree->GetEntry(0);
  REQUIRE(readcdata->Data.size() == 1);

  tree->GetEntry(1);
  REQUIRE(readcdata->Data.size() == 2);

}
