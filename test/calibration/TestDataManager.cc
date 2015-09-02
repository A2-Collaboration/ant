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
void dotest_extendable();

TEST_CASE("CalibrationDataManager: Save/Load","[calibration]")
{
    tmpfolder_t tmp;
    auto ndata = dotest_store(tmp.foldername);
    dotest_load(tmp.foldername,ndata);
    dotest_changes(tmp.foldername);
}

TEST_CASE("CalibrationDataManager: Extendable","[calibration]")
{
    dotest_extendable();
}

void dotest_extendable()
{
    // some commonly used variables
    tmpfolder_t tmp;
    TCalibrationData cdata;
    string id;
    DataManager calibman(tmp.foldername);
    auto check = [] (const TCalibrationData& c, const interval<TID>& i) {
        return c.FirstID == i.Start() && c.LastID == i.Stop();
    };
    const interval<TID> i1(TID(4),  TID(6));
    const interval<TID> i2(TID(9),  TID(10));
    const interval<TID> i3(TID(12), TID(15));
    const interval<TID> i4(TID(2),  TID(16));

    // first simple test
    id = "1";
    calibman.Add(TCalibrationData(id, i1.Start(), i1.Stop(), true));
    calibman.Add(TCalibrationData(id, i2.Start(), i2.Stop(), true));
    calibman.Add(TCalibrationData(id, i3.Start(), i3.Stop(), true));

    REQUIRE(calibman.GetData(id, TID(3), cdata));
    REQUIRE(check(cdata, i1));

    REQUIRE(calibman.GetData(id, TID(7), cdata));
    REQUIRE(check(cdata, i2));

    REQUIRE(calibman.GetData(id, TID(11), cdata));
    REQUIRE(check(cdata, i3));

    REQUIRE(calibman.GetData(id, TID(16), cdata));
    REQUIRE(check(cdata, i3));


    id = "2";
    calibman.Add(TCalibrationData(id, i4.Start(), i4.Stop(), true));
    calibman.Add(TCalibrationData(id, i1.Start(), i1.Stop(), true));
    calibman.Add(TCalibrationData(id, i2.Start(), i2.Stop(), true));
    calibman.Add(TCalibrationData(id, i3.Start(), i3.Stop(), true));

    REQUIRE(calibman.GetData(id, TID(3), cdata));
    REQUIRE(check(cdata, i4));

    REQUIRE(calibman.GetData(id, TID(1), cdata));
    REQUIRE(check(cdata, i4));

    REQUIRE(calibman.GetData(id, TID(17), cdata));
    REQUIRE(check(cdata, i4));

    REQUIRE(calibman.GetData(id, TID(9), cdata));
    REQUIRE(check(cdata, i2));
}

unsigned dotest_store(const string& foldername)
{


    DataManager calibman(foldername);

    TCalibrationData cdata("M",
                           "comment",
                           0,
                           "1",
                           TID(0),TID(16)
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
                             TID(first),TID(last)
                             );
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

    calibman.GetData("1",TID(0),cdata);
    REQUIRE(cdata.TimeStamp == 0);

    calibman.GetData("1",TID(1),cdata);
    REQUIRE(cdata.TimeStamp == 0);

    calibman.GetData("1",TID(3),cdata);
    REQUIRE(cdata.TimeStamp == 3);

    calibman.GetData("1",TID(4),cdata);
    REQUIRE(cdata.TimeStamp == 3);

    calibman.GetData("1",TID(5),cdata);
    REQUIRE(cdata.TimeStamp == 4);

    calibman.GetData("1",TID(14),cdata);
    REQUIRE(cdata.TimeStamp == 7);

    REQUIRE_FALSE(calibman.GetData("1",TID(21),cdata));

    calibman.GetData("1",TID(23),cdata);
    REQUIRE(cdata.TimeStamp == 6);

    REQUIRE_FALSE(calibman.GetData("1",TID(29),cdata));

    auto genchange = calibman.GetChangePoints("1");

    list<TID> manchange({TID(0),
                         TID(2),
                         TID(3),
                         TID(5),
                         TID(8),
                         TID(9),
                         TID(13),
                         TID(14),
                         TID(15),
                         TID(21),
                         TID(22),
                         TID(25)
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
