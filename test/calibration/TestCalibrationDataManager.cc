#include "catch.hpp"

#include "tmpfile_t.h"
#include "DataManager.h"
#include "TCalibrationData.h"
#include "TDataRecord.h"

#include <list>
#include <algorithm>


using namespace std;
using namespace ant;
using namespace ant::calibration;

unsigned dotest_store(const string& filename);
void dotest_load(const string& filename, unsigned ndata);
void dotest_changes(const string& filename);

TEST_CASE("CalibrationDataManager","[calibration]")
{
    tmpfile_t tmpfile;
    auto ndata = dotest_store(tmpfile.filename);
    dotest_load(tmpfile.filename,ndata);
    dotest_changes(tmpfile.filename);
}

unsigned dotest_store(const string& filename)
{


    DataManager calibman(filename);

    TCalibrationData cdata("M",
                           "comment",
                           0,
                           "1",
                           TID(0,0),TID(0,16)
                           );
    cdata.Data.emplace_back(0,1);
    cdata.Data.emplace_back(1,2);
    calibman.Add(cdata);
    unsigned ndata(1);

    auto mdata = [&cdata,&ndata] (unsigned first, unsigned last, unsigned time)
    {
        ndata++;
        TCalibrationData tmp(cdata.Author,
                                cdata.Comment,
                                time,
                                cdata.CalibrationID,
                                TID(0,first),TID(0,last));
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
    REQUIRE(calibman.GetNumberOfDataPoints("1") == ndata);
    REQUIRE(calibman.GetNumberOfDataPoints("2") == 3);

    return ndata;
}

void dotest_load(const string &filename,unsigned ndata)
{
    DataManager calibman(filename);
    REQUIRE(calibman.GetNumberOfCalibrations() == 2);
    REQUIRE(calibman.GetNumberOfDataPoints("1") == ndata);
    REQUIRE(calibman.GetNumberOfDataPoints("2") == 3);
}

void dotest_changes(const string& filename)
{
    DataManager calibman(filename);

    TCalibrationData cdata;
//    interval<TID> idRangeTEST(TID(0,0),TID(0,24));
//    interval<TID> idRange(TID(1,1),TID(1,1));

//    REQUIRE(calibman.GetIDRange("1",idRange));
//    REQUIRE( idRangeTEST == idRange);

    calibman.GetData("1",TID(0,0),cdata);
    REQUIRE(cdata.TimeStamp == 0);

    calibman.GetData("1",TID(0,1),cdata);
    REQUIRE(cdata.TimeStamp == 0);

    calibman.GetData("1",TID(0,3),cdata);
    REQUIRE(cdata.TimeStamp == 3);

    calibman.GetData("1",TID(0,4),cdata);
    REQUIRE(cdata.TimeStamp == 3);

    calibman.GetData("1",TID(0,5),cdata);
    REQUIRE(cdata.TimeStamp == 4);

    calibman.GetData("1",TID(0,14),cdata);
    REQUIRE(cdata.TimeStamp == 7);

    REQUIRE_FALSE(calibman.GetData("1",TID(0,21),cdata));

    calibman.GetData("1",TID(0,23),cdata);
    REQUIRE(cdata.TimeStamp == 6);

    REQUIRE_FALSE(calibman.GetData("1",TID(0,29),cdata));

    auto genchange = calibman.GetChangePoints("1");

    list<TID> manchange({TID(0,0),
                           TID(0,2),
                           TID(0,3),
                           TID(0,5),
                           TID(0,8),
                           TID(0,9),
                           TID(0,13),
                           TID(0,14),
                           TID(0,15),
                           TID(0,21),
                           TID(0,22),
                           TID(0,25)
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
