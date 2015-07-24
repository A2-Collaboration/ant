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

    calibman.Add(TCalibrationData("M",
                                  "comment",
                                  1234567890,
                                  "setupID",
                                  TID(1,2),TID(1,3),
                                  {{0,1},{1,2}}));
}

void dotest_load(const string &filename)
{
    REQUIRE_NOTHROW(
                CalibrationManager calibman(filename)
                );
}
