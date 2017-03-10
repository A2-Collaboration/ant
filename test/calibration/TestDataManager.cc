#include "catch.hpp"

#include "DataManager.h"
#include "DataBase.h"

#include "tree/TCalibrationData.h"

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
    calibman.Add(cdata,  Calibration::AddMode_t::AsDefault);
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

    calibman.Add(mdata( 4,  4, 1), Calibration::AddMode_t::StrictRange);
    REQUIRE_THROWS_AS(calibman.Add(mdata( 2,  108, 2), Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    REQUIRE_THROWS_AS(calibman.Add(mdata( 3,  106, 3), Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    calibman.Add(mdata( 5,  7, 4), Calibration::AddMode_t::StrictRange);
    calibman.Add(mdata(13, 20, 5), Calibration::AddMode_t::StrictRange);
    calibman.Add(mdata(22, 24, 6), Calibration::AddMode_t::StrictRange);
    REQUIRE_THROWS_AS(calibman.Add(mdata(14, 124, 7), Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);

    // three times the Add failed above
    ndata -= 3;

    // make some simpler second calibration without defaults
    cdata.CalibrationID = "2";
    cdata.TimeStamp = 0;
    cdata.FirstID.Lower = 1;
    cdata.LastID.Lower = 1;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 3;
    cdata.LastID.Lower = 6;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 7;
    cdata.LastID.Lower = 7;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    // make something with more complicated time ranges
    cdata.CalibrationID = "3";
    cdata.TimeStamp = 0;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Lower = 8;
    calibman.Add(cdata, Calibration::AddMode_t::AsDefault);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Lower = 0xffffffff;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 100000;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Timestamp = cdata.FirstID.Timestamp;
    cdata.LastID.Lower = 3;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    // add things multiple times
    cdata.CalibrationID = "4";
    cdata.TimeStamp = 0;
    cdata.FirstID.Timestamp = 0;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Timestamp = cdata.FirstID.Timestamp;
    cdata.LastID.Lower = 8;
    calibman.Add(cdata, Calibration::AddMode_t::AsDefault);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::AsDefault);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::AsDefault);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::AsDefault);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::AsDefault);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Lower = 1;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 200000;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Timestamp = cdata.FirstID.Timestamp;
    cdata.LastID.Lower = 3;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);
    cdata.TimeStamp++;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    // test RightOpen intervals
    cdata.CalibrationID = "5";
    cdata.TimeStamp = 0;
    cdata.FirstID.Timestamp = 0;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Timestamp = cdata.FirstID.Timestamp;
    cdata.LastID.Lower = 8;
    calibman.Add(cdata, Calibration::AddMode_t::AsDefault);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 10;
    cdata.FirstID.Lower = 0;
    calibman.Add(cdata, Calibration::AddMode_t::RightOpen);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 20;
    cdata.FirstID.Lower = 10;
    calibman.Add(cdata, Calibration::AddMode_t::RightOpen);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 5;
    cdata.FirstID.Lower = 0;
    calibman.Add(cdata, Calibration::AddMode_t::RightOpen);


    // test RightOpen intervals more complicated
    cdata.CalibrationID = "6";
    cdata.TimeStamp = 0;
    cdata.FirstID.Timestamp = 10;
    cdata.FirstID.Lower = 0;
    calibman.Add(cdata, Calibration::AddMode_t::RightOpen);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 20;
    cdata.FirstID.Lower = 10;
    calibman.Add(cdata, Calibration::AddMode_t::RightOpen);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 5;
    cdata.FirstID.Lower = 0;
    calibman.Add(cdata, Calibration::AddMode_t::RightOpen);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 10;
    cdata.FirstID.Lower = 0;
    calibman.Add(cdata, Calibration::AddMode_t::RightOpen);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 5;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Timestamp = 9;
    cdata.LastID.Lower = 0xffffffff;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    // test AdHoc Data intervals (cannot be loaded from database, issues warning)
    // but are still stored there
    cdata.CalibrationID = "7";
    cdata.TimeStamp = 0;
    cdata.FirstID = TID(0, 0, {TID::Flags_t::AdHoc});
    cdata.LastID = TID(0, 5, {TID::Flags_t::AdHoc});
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 5;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Timestamp = 9;
    cdata.LastID.Lower = 0xffffffff;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 5;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Timestamp = 10;
    cdata.LastID.Lower = 0xffffffff;
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);


    // test AdHoc MC intervals (should always be added to default)
    cdata.CalibrationID = "8";
    cdata.TimeStamp = 0;
    cdata.FirstID = TID(0, 0, {TID::Flags_t::MC});
    cdata.LastID = TID(0, 5, {TID::Flags_t::MC});
    calibman.Add(cdata, Calibration::AddMode_t::RightOpen);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 5;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Timestamp = 9;
    cdata.LastID.Lower = 0xffffffff;
    calibman.Add(cdata, Calibration::AddMode_t::AsDefault);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 5;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Timestamp = 9;
    cdata.LastID.Lower = 0xffffffff;
    calibman.Add(cdata, Calibration::AddMode_t::StrictRange);

    cdata.TimeStamp++;
    cdata.FirstID.Timestamp = 5;
    cdata.FirstID.Lower = 0;
    cdata.LastID.Timestamp = 11;
    cdata.LastID.Lower = 0xffffffff;
    REQUIRE_NOTHROW(calibman.Add(cdata, Calibration::AddMode_t::StrictRange));

    // some not allowed things
    cdata.FirstID = TID();
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::RightOpen),
                      DataBase::Exception);
    cdata.FirstID = TID(0,0u);
    cdata.LastID = TID();
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    cdata.FirstID = TID(0,10u);
    cdata.LastID = TID(0,0u);
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    cdata.FirstID = TID(0,0u,{TID::Flags_t::MC});
    cdata.LastID = TID(0,10u);
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    cdata.FirstID = TID(0,0u,{TID::Flags_t::AdHoc});
    cdata.LastID = TID(0,10u);
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    cdata.FirstID = TID(0,0u);
    cdata.LastID = TID(0,10u,{TID::Flags_t::MC});
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    cdata.FirstID = TID(0,0u);
    cdata.LastID = TID(0,10u,{TID::Flags_t::AdHoc});
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);
    cdata.CalibrationID = "";
    REQUIRE_THROWS_AS(calibman.Add(cdata, Calibration::AddMode_t::StrictRange),
                      DataBase::Exception);


    // check status
    REQUIRE(calibman.GetNumberOfCalibrationIDs() == 8);
    REQUIRE(calibman.GetNumberOfCalibrationData("1") == ndata);
    REQUIRE(calibman.GetNumberOfCalibrationData("2") == 3);
    REQUIRE(calibman.GetNumberOfCalibrationData("3") == 3);
    REQUIRE(calibman.GetNumberOfCalibrationData("4") == 12);
    REQUIRE(calibman.GetNumberOfCalibrationData("5") == 4);
    REQUIRE(calibman.GetNumberOfCalibrationData("6") == 5);
    REQUIRE(calibman.GetNumberOfCalibrationData("7") == 2); // one failed intentionally
    REQUIRE(calibman.GetNumberOfCalibrationData("8") == 4);



    return ndata;
}

void dotest_load(const string &foldername,unsigned ndata)
{
    DataManager calibman(foldername);
    REQUIRE(calibman.GetNumberOfCalibrationIDs() == 8);
    REQUIRE(calibman.GetNumberOfCalibrationData("1") == ndata);
    REQUIRE(calibman.GetNumberOfCalibrationData("2") == 3);
    REQUIRE(calibman.GetNumberOfCalibrationData("3") == 3);
    REQUIRE(calibman.GetNumberOfCalibrationData("4") == 12);
    REQUIRE(calibman.GetNumberOfCalibrationData("5") == 4);
    REQUIRE(calibman.GetNumberOfCalibrationData("6") == 5);
    REQUIRE(calibman.GetNumberOfCalibrationData("7") == 2);
    REQUIRE(calibman.GetNumberOfCalibrationData("8") == 4);
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

    // test 4
    REQUIRE(calibman.GetData("4",TID(0,2u),cdata));
    REQUIRE(cdata.TimeStamp == 4);

    REQUIRE(calibman.GetData("4",TID(0,0u),cdata));
    REQUIRE(cdata.TimeStamp == 7);

    REQUIRE(calibman.GetData("4",TID(200000,2u),cdata));
    REQUIRE(cdata.TimeStamp == 11);

    // test 5
    REQUIRE(calibman.GetData("5",TID(0,2u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 0);
    REQUIRE(nextChangePoint == TID(5,0u));

    REQUIRE(calibman.GetData("5",TID(5,0u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 3);
    REQUIRE(nextChangePoint == TID(10,0u));

    REQUIRE(calibman.GetData("5",TID(10,0u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 1);
    REQUIRE(nextChangePoint == TID(20,10u));

    REQUIRE(calibman.GetData("5",TID(20,10u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 0);
    REQUIRE(nextChangePoint.IsInvalid());

    // test 6
    REQUIRE_FALSE(calibman.GetData("6",TID(0,2u),cdata,nextChangePoint));

    REQUIRE(calibman.GetData("6",TID(5,10u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 4);
    REQUIRE(nextChangePoint == TID(10,0u));

    REQUIRE(calibman.GetData("6",TID(10,0u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 3);
    REQUIRE(nextChangePoint == TID(20,10u));

    REQUIRE(calibman.GetData("6",TID(20,20u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 1);
    REQUIRE(nextChangePoint.IsInvalid());

    // test 7
    // Data TID with AdHoc flag don't load anything
    REQUIRE_FALSE(calibman.GetData("7",TID(0,2u,{TID::Flags_t::AdHoc}),cdata,nextChangePoint));
    REQUIRE_FALSE(calibman.GetData("7",TID(40,0u),cdata,nextChangePoint));
    // but omitting the AdHoc flag shows the data
    REQUIRE(calibman.GetData("7",TID(6,0u),cdata,nextChangePoint));
    REQUIRE(cdata.TimeStamp == 1);

    // test 8
    // MC flag always gets last added
    REQUIRE(calibman.GetData("8",TID(0,2u,{TID::Flags_t::MC}),cdata,nextChangePoint));
    REQUIRE(nextChangePoint.IsInvalid());

    REQUIRE(cdata.FirstID.Timestamp == 5);
    REQUIRE(cdata.FirstID.Lower == 0);
    REQUIRE(cdata.LastID.Timestamp == 11);
    REQUIRE(cdata.LastID.Lower == 0xffffffff);


}
