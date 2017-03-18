#include "catch.hpp"

#include "tmpfile_t.h"

#include "tree/TCalibrationData.h"

#include "base/WrapTTree.h"
#include "base/WrapTFile.h"

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("TCalibrationData: Write to TTree", "[tree]") {
    dotest();
}

void dotest()
{
    tmpfile_t tmpfile;

    struct Tree_t : WrapTTree {
        ADD_BRANCH_T(TCalibrationData, cdata)
        ADD_BRANCH_T(TCalibrationData, cdata2)
    };

    {
        WrapTFileOutput outputfile(tmpfile.filename);
        Tree_t t;
        REQUIRE_NOTHROW(t.CreateBranches(outputfile.CreateInside<TTree>("test","")));
        t.cdata = TCalibrationData(tmpfile.filename,
                                   ant::TID(0,0u),ant::TID(0,1u,{ant::TID::Flags_t::MC}));

        t.cdata2 = TCalibrationData(tmpfile.filename,
                                    ant::TID(0,1011u),
                                    ant::TID(0,1213u,{ant::TID::Flags_t::MC})
                                    );

        t.cdata2().Data.emplace_back(1,1.);
        t.cdata2().Data.emplace_back(2,2.1);

        t.cdata().Author  = "Martin Wolfes";


        t.cdata().Data.emplace_back(1,2.1);
        t.Tree->Fill();
        t.cdata().Data.emplace_back(2,3.2);
        t.Tree->Fill();
    }

    {
        WrapTFileInput inputfile(tmpfile.filename);

        Tree_t t;
        REQUIRE(inputfile.GetObject("test",t.Tree));
        REQUIRE_NOTHROW(t.LinkBranches());

        REQUIRE(t.Tree->GetEntries() == 2);

        t.Tree->GetEntry(0);
        REQUIRE(t.cdata().Data.size() == 1);

        t.Tree->GetEntry(1);
        REQUIRE(t.cdata().Data.size() == 2);
        REQUIRE(t.cdata().Data[0].Key == 1);
        REQUIRE(t.cdata().Data[1].Value == 3.2);

    }

}
