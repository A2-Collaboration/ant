#include "catch.hpp"

#include "base/WrapTTree.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"
#include "base/std_ext/memory.h"

#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;
using namespace ant;

void dotest();
void dotest_copy();


TEST_CASE("WrapTTree: Basics", "[base]") {
    dotest();
}

TEST_CASE("WrapTTree: Copy", "[base]") {
    dotest_copy();
}

struct MyTree : WrapTTree {
    ADD_BRANCH_T(bool,           Flag1)       // simple type
    ADD_BRANCH_T(unsigned,       N1)          // simple type
    ADD_BRANCH_T(short,          N2)          // simple type
    ADD_BRANCH_T(bool,           Flag2, true) // simple type
    ADD_BRANCH_T(vector<double>, Array, 3)    // use ctor of std::vector
    ADD_BRANCH_T(TLorentzVector, LV)          // complex ROOT type
};

void dotest() {

    tmpfile_t tmpfile;

    // do some simple IO
    {
        WrapTFileOutput outputfile(tmpfile.filename);
        outputfile.cd();

        MyTree t;
        t.CreateBranches(outputfile.CreateInside<TTree>("test","test"));
        REQUIRE(t.Tree != nullptr);

        REQUIRE(t.Flag2);
        t.N1 = 3;
        t.N2 = 5;
        REQUIRE(t.Array().size() == 3);
        t.LV = TLorentzVector(1,2,3,4);
        t.Tree->Fill();
        REQUIRE(t.Tree->GetEntries()==1);
    }
    {
        WrapTFileInput inputfile(tmpfile.filename);
        MyTree t;
        REQUIRE(inputfile.GetObject("test",t.Tree));
        t.LinkBranches(t.Tree);
        REQUIRE(t.Tree->GetEntries()==1);
        t.Tree->GetEntry(0);
        REQUIRE(t.N1 == 3);
        REQUIRE(t.N2 == 5);
        REQUIRE(t.Array().size() == 3);
        REQUIRE(t.LV().E() == Approx(4));
    }

}

void dotest_copy() {

    MyTree t1;
    t1.N1 = 5;
    t1.N2 = 7;
    t1.Array().at(2) = 3;
    t1.LV = TLorentzVector(1,2,3,4);

    // test with same type
    {
        MyTree t2;
        t2 = t1;
        REQUIRE(t2.N1 == 5);
        REQUIRE(t2.N2 == 7);
        REQUIRE(t2.Array().at(2) == 3);
        REQUIRE(t2.LV().E() == Approx(4));
    }

    // test with derived type
    {
        REQUIRE_FALSE(t1.Flag1);
        t1.Flag1 = true;
        struct MyTree2 : MyTree {
            ADD_BRANCH_T(double, SomeValue)
        };

        MyTree2 t2;
        REQUIRE(t2.CopyFrom(t1));

        REQUIRE(t2.Flag1);
        REQUIRE(t2.Array().at(2) == 3);
        REQUIRE(t2.LV().E() == Approx(4));
    }
}
