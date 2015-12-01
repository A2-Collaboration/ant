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
    calibman.Add(cdata);
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

    calibman.Add(mdata( 4,  4, 1));
    calibman.Add(mdata( 2,  8, 2));
    calibman.Add(mdata( 3,  6, 3));
    calibman.Add(mdata( 5,  7, 4));
    calibman.Add(mdata(13, 20, 5));
    calibman.Add(mdata(22, 24, 6));
    calibman.Add(mdata(14, 14, 7));

    cdata.CalibrationID = "2";
    cdata.TimeStamp++;
    cdata.FirstID.Lower = 2;
    cdata.LastID.Lower = 8;
    calibman.Add(cdata);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 3;
    cdata.LastID.Lower = 6;
    calibman.Add(cdata);

    cdata.TimeStamp++;
    cdata.FirstID.Lower = 5;
    cdata.LastID.Lower = 7;
    calibman.Add(cdata);

    REQUIRE(calibman.GetNumberOfCalibrations() == 2);
    REQUIRE(calibman.GetNumberOfDataPoints("1") == ndata);
    REQUIRE(calibman.GetNumberOfDataPoints("2") == 3);

    return ndata;
}

void dotest_load(const string &foldername,unsigned ndata)
{
    DataManager calibman(foldername);
    REQUIRE(calibman.GetNumberOfCalibrations() == 2);
    REQUIRE(calibman.GetNumberOfDataPoints("1") == ndata);
    REQUIRE(calibman.GetNumberOfDataPoints("2") == 3);
}

void dotest_changes(const string& foldername)
{
    DataManager calibman(foldername);

    TCalibrationData cdata;
    //    interval<TID> idRangeTEST(TID(0,0),TID(0,24));
    //    interval<TID> idRange(TID(1,1),TID(1,1));

    //    REQUIRE(calibman.GetIDRange("1",idRange));
    //    REQUIRE( idRangeTEST == idRange);

    calibman.GetData("1",TID(0,0u),cdata);
    REQUIRE(cdata.TimeStamp == 0);

    calibman.GetData("1",TID(0,1u),cdata);
    REQUIRE(cdata.TimeStamp == 0);

    calibman.GetData("1",TID(0,3u),cdata);
    REQUIRE(cdata.TimeStamp == 3);

    calibman.GetData("1",TID(0,4u),cdata);
    REQUIRE(cdata.TimeStamp == 3);

    calibman.GetData("1",TID(0,5u),cdata);
    REQUIRE(cdata.TimeStamp == 4);

    calibman.GetData("1",TID(0,14u),cdata);
    REQUIRE(cdata.TimeStamp == 7);

    REQUIRE_FALSE(calibman.GetData("1",TID(0,21u),cdata));

    calibman.GetData("1",TID(0,23u),cdata);
    REQUIRE(cdata.TimeStamp == 6);

    REQUIRE_FALSE(calibman.GetData("1",TID(0,29u),cdata));

    auto genchange = calibman.GetChangePoints("1");

    list<TID> manchange({TID(0,0u),
                         TID(0,2u),
                         TID(0,3u),
                         TID(0,5u),
                         TID(0,8u),
                         TID(0,9u),
                         TID(0,13u),
                         TID(0,14u),
                         TID(0,15u),
                         TID(0,21u),
                         TID(0,22u),
                         TID(0,25u)
                        });

    REQUIRE(genchange.size() == manchange.size());

    auto itg = genchange.begin();
    auto itm = manchange.begin();

    while (itg != genchange.end())
    {
        REQUIRE(*itg == *itm);
        itg++;
        itm++;
    }

}
