#include "catch.hpp"

#include "tmpfile_t.h"
#include "CalibrationManager.h"
#include "TCalibrationData.h"
#include "TDataRecord.h"

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("CalibrationManager","[calibration]")
{
    dotest();
}


void dotest()
{
    ant::tmpfile_t tmpfile;

    CalibrationManager calibman(tmpfile.filename);

    REQUIRE_NOTHROW(calibman.Add(TCalibrationData("M",
                                  "comment",
                                  1234567890,
                                  "setupID",
                                  TID(1,2),TID(1,3),
                                  {{0,1},{1,2}})));








}
