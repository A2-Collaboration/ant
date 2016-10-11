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
    ADD_BRANCH_T(unsigned,       A)    // simple type
    ADD_BRANCH_T(vector<double>, B, 3) // use ctor of std::vector
    ADD_BRANCH_T(TLorentzVector, C)    // complex ROOT type
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

        t.A = 3;
        REQUIRE(t.B().size() == 3);
        t.C = TLorentzVector(1,2,3,4);
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
        REQUIRE(t.B().size() == 3);
        REQUIRE(t.C().E() == Approx(4));
    }

}

void dotest_copy() {

    MyTree t1;
    t1.A = 5;
    t1.B().at(2) = 3;
    t1.C = TLorentzVector(1,2,3,4);

    // test with same type
    {
        MyTree t2;
        t2 = t1;
        REQUIRE(t2.A == 5);
        REQUIRE(t2.C().E() == Approx(4));
    }

//    struct MyTree2 : MyTree {
//        ADD_BRANCH_T(double, D)
//    };

//    MyTree2 t3;
//    t3 = t2;
}
