#include "catch.hpp"

#include "DataManager.h"

#include "tree/TCalibrationData.h"
#include "tree/TDataRecord.h"

#include "base/tmpfile_t.h"
#include "base/interval.h"

#include <list>
#include <algorithm>


using namespace std;
using namespace ant;
using namespace ant::calibration;

unsigned dotest_store(const string& foldername);
void dotest_load(const string& foldername, unsigned ndata);
void dotest_changes(const string& foldername);

TEST_CASE("CalibrationDataManager: Save/Load","[calibration]")
{
    tmpfolder_t tmp;
    auto ndata = dotest_store(tmp.foldername);
    dotest_load(tmp.foldername,ndata);
    dotest_changes(tmp.foldername);
}

unsigned dotest_store(const string& foldername)
{
    DataManager calibman(foldername);

    TCalibrationData cdata("1",
                           TID(0,0u),TID(0,16u)
                           );
    cdata.TimeStamp = 0;
    cdata.Data.emplace_back(0,1);
    cdata.Data.emplace_back(1,2);
    calibman.Add(cdata,  DataBase::mode_t::AsDefault);
    unsigned ndata(1);

    auto mdata = [&cdata,&ndata] (unsigned first, unsigned last, unsigned time)
    {
        ndata++;
        TCalibrationData tmp(cdata.CalibrationID,
                             TID(0,first),TID(0,last)
                             );
        tmp.Author = cdata.Author;
        tmp.TimeStamp = time;
        tmp.Data = cdata.Data;
        return tmp;
    };

    calibman.Add(mdata( 4,  4, 1), DataBase::mode_t::StrictRange);
    REQUIRE_THROWS_AS(calibman.Add(mdata( 2,  8, 2), DataBase::mode_t::StrictRange), DataBase::Exception);
    REQUIRE_THROWS_AS(calibman.Add(mdata( 3,  6, 3), DataBase::mode_t::StrictRange), DataBase::Exception);
    calibman.Add(mdata( 5,  7, 4), DataBase::mode_t::StrictRange);
    calibman.Add(mdata(13, 20, 5), DataBase::mode_t::StrictRange);
    calibman.Add(mdata(22, 24, 6), DataBase::mode_t::StrictRange);
    REQUIRE_THROWS_AS(calibman.Add(mdata(14, 14, 7), DataBase::mode_t::StrictRange), DataBase::Exception);

    // three times the Add failed above
    ndata -= 3;

    // make some simpler second calibration without defaults
    cdata.CalibrationID = "2";
    cdata.TimeStamp = 0;
    cdata.FirstID.Lower = 1;
    cdata.LastID.Lower = 1;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 3;
    cdata.LastID.Lower = 6;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 7;
    cdata.LastID.Lower = 7;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);

    // make something with more complicated time ranges
    cdata.CalibrationID = "3";
    cdata.TimeStamp = 0;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Lower = 8;
    calibman.Add(cdata, DataBase::mode_t::AsDefault);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Lower = 0xffffffff;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 100000;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Timestamp = cdata.FirstID.Timestamp;
    cdata.LastID.Lower = 3;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);

    // add things multiple times
    cdata.CalibrationID = "4";
    cdata.TimeStamp = 0;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Lower = 8;
    calibman.Add(cdata, DataBase::mode_t::AsDefault);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::AsDefault);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::AsDefault);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::AsDefault);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::AsDefault);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Lower = 1;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 12302193;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Timestamp = cdata.FirstID.Timestamp;
    cdata.LastID.Lower = 3;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, DataBase::mode_t::StrictRange);

    REQUIRE(calibman.GetNumberOfCalibrationIDs() == 4);
    REQUIRE(calibman.GetNumberOfCalibrationData("1") == ndata);
    REQUIRE(calibman.GetNumberOfCalibrationData("2") == 3);
    REQUIRE(calibman.GetNumberOfCalibrationData("3") == 3);
    REQUIRE(calibman.GetNumberOfCalibrationData("4") == 12);

    return ndata;
}

void dotest_load(const string &foldername,unsigned ndata)
{
    DataManager calibman(foldername);
    REQUIRE(calibman.GetNumberOfCalibrationIDs() == 4);
    REQUIRE(calibman.GetNumberOfCalibrationData("1") == ndata);
    REQUIRE(calibman.GetNumberOfCalibrationData("2") == 3);
    REQUIRE(calibman.GetNumberOfCalibrationData("3") == 3);
    REQUIRE(calibman.GetNumberOfCalibrationData("4") == 12);
}

void dotest_changes(const string& foldername)
{
    DataManager calibman(foldername);

    TCalibrationData cdata;
    TID nextChangePoint;

    // test 1

    REQUIRE(calibman.GetData("1", TID(0,0u), cdata));
    REQUIRE(cdata.TimeStamp == 0);

    REQUIRE(calibman.GetData("1",TID(0,1u),cdata));
    REQUIRE(cdata.TimeStamp == 0);

    REQUIRE(calibman.GetData("1",TID(0,3u),cdata));
    REQUIRE(cdata.TimeStamp == 0);

    REQUIRE(calibman.GetData("1",TID(0,4u),cdata));
    REQUIRE(cdata.TimeStamp == 1);

    REQUIRE(calibman.GetData("1",TID(0,5u),cdata));
    REQUIRE(cdata.TimeStamp == 4);

    REQUIRE(calibman.GetData("1",TID(0,14u),cdata));
    REQUIRE(cdata.TimeStamp == 5);

    REQUIRE(calibman.GetData("1",TID(0,21u),cdata));
    REQUIRE(cdata.TimeStamp == 0);

    REQUIRE(calibman.GetData("1",TID(0,23u),cdata));
    REQUIRE(cdata.TimeStamp == 6);

    REQUIRE(calibman.GetData("1",TID(0,26u),cdata));
    REQUIRE(cdata.TimeStamp == 0);

    REQUIRE(calibman.GetData("1",TID(0,0u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 0);
    REQUIRE(nextChangePoint == TID(0,4u));

    // test 2
    REQUIRE_FALSE(calibman.GetData("2",TID(0,0u),cdata,nextChangePoint));
    REQUIRE(nextChangePoint == TID(0,1u));

    REQUIRE(calibman.GetData("2",TID(0,7u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 2);
    REQUIRE_FALSE(nextChangePoint.IsInvalid());
    REQUIRE(nextChangePoint == TID(0,8u));

    REQUIRE_FALSE(calibman.GetData("2",TID(0,10u),cdata,nextChangePoint));
    REQUIRE(nextChangePoint.IsInvalid());

     // test 3
    REQUIRE(calibman.GetData("3",TID(0,0u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 1);
    REQUIRE(nextChangePoint == TID(1, 0u)); // check "overflow" of TID

    REQUIRE(calibman.GetData("3",TID(100000,2u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 2);
}
