#include "catch.hpp"

#include "tmpfile_t.h"
#include "CalibrationManager.h"
#include "TCalibrationData.h"
#include "TDataRecord.h"

using namespace std;
using namespace ant;

void dotest_store(const string& filename);
void dotest_load(const string& filename);

TEST_CASE("CalibrationManager","[calibration]")
{
    tmpfile_t tmpfile;
    dotest_store(tmpfile.filename);
    dotest_load(tmpfile.filename);
}


void dotest_store(const string& filename)
{

    CalibrationManager calibman(filename);

    TCalibrationData cdata("M",
                           "comment",
                           std::time(nullptr),
                           "1",
                           TID(0,0),TID(0,10),
                           {{0,1},{1,2}});

    calibman.Add(cdata);

    cdata.TimeStamp++;
    cdata.FirstID.Value = 2;
    cdata.LastID.Value = 8;
    calibman.Add(cdata);

    cdata.TimeStamp++;
    cdata.FirstID.Value = 3;
    cdata.LastID.Value = 6;
    calibman.Add(cdata);

    cdata.TimeStamp++;
    cdata.FirstID.Value = 5;
    cdata.LastID.Value = 7;
    calibman.Add(cdata);

    cdata.SetupID = "2";
    cdata.TimeStamp++;
    cdata.FirstID.Value = 2;
    cdata.LastID.Value = 8;
    calibman.Add(cdata);

    cdata.TimeStamp++;
    cdata.FirstID.Value = 3;
    cdata.LastID.Value = 6;
    calibman.Add(cdata);

    cdata.TimeStamp++;
    cdata.FirstID.Value = 5;
    cdata.LastID.Value = 7;
    calibman.Add(cdata);

    REQUIRE(calibman.GetNumberOfCalibrations() == 2);
    REQUIRE(calibman.GetNumberOfDataPoints("1") == 4);
    REQUIRE(calibman.GetNumberOfDataPoints("2") == 3);


}

void dotest_load(const string &filename)
{
    CalibrationManager calibman(filename);
    REQUIRE(calibman.GetNumberOfCalibrations() == 2);
    REQUIRE(calibman.GetNumberOfDataPoints("1") == 4);
    REQUIRE(calibman.GetNumberOfDataPoints("2") == 3);
}
