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
void dotest_nasty();
void dotest_pod();


TEST_CASE("WrapTTree: Basics", "[base]") {
    dotest();
}

TEST_CASE("WrapTTree: Copy", "[base]") {
    dotest_copy();
}

TEST_CASE("WrapTTree: Nasty", "[base]") {
    dotest_nasty();
}

TEST_CASE("WrapTTree: Plain-old datatypes", "[base]") {
    dotest_pod();
}


struct MyTree : WrapTTree {
    ADD_BRANCH_T(bool,           Flag1)        // simple type
    ADD_BRANCH_T(unsigned,       N1)           // simple type
    ADD_BRANCH_T(short,          N2)           // simple type
    ADD_BRANCH_T(bool,           Flag2, true)  // simple type
    ADD_BRANCH_T(vector<double>, SomeArray, 3) // use ctor of std::vector
    ADD_BRANCH_T(TLorentzVector, LV)           // complex ROOT type
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
        REQUIRE(t.SomeArray().size() == 3);
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
        REQUIRE(t.SomeArray().size() == 3);
        REQUIRE(t.LV().E() == Approx(4));
    }

}

void dotest_copy() {

    MyTree t1;
    t1.N1 = 5;
    t1.N2 = 7;
    t1.SomeArray().at(2) = 3;
    t1.LV = TLorentzVector(1,2,3,4);

    // test with same type
    {
        MyTree t2;
        t2 = t1;
        REQUIRE(t2.N1 == 5);
        REQUIRE(t2.N2 == 7);
        REQUIRE(t2.SomeArray().at(2) == 3);
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
        REQUIRE(t2.SomeArray().at(2) == 3);
        REQUIRE(t2.LV().E() == Approx(4));
    }
}

void dotest_nasty() {
    {
        // same name
        struct MyTree2 : MyTree {
            ADD_BRANCH_T(double, SomeArray)
        };
        REQUIRE_THROWS_AS(std_ext::make_unique<MyTree2>(),WrapTTree::Exception);
    }

    {
        // nullptr to create branches
        MyTree t;
        REQUIRE_THROWS_AS(t.CreateBranches(nullptr), WrapTTree::Exception);
    }
}

void dotest_pod() {
    constexpr auto nFills = 100; // not higher than max, see below

    tmpfile_t tmpfile;
    std::vector<double> filled_data;
    {
        WrapTFileOutput outputfile(tmpfile.filename, WrapTFileOutput::mode_t::recreate, true);
        // the following code is copied from
        // https://github.com/A2-Collaboration-dev/acqu/blob/master/acqu_user/root/src/TA2GoAT.h
        // https://github.com/A2-Collaboration-dev/acqu/blob/master/acqu_user/root/src/TA2GoAT.cc
        auto treeTracks = new TTree("treeTracks","treeTracks");
        constexpr auto TA2GoAT_MAX_PARTICLE = 128; // was #define TA2GoAT_MAX_PARTICLE 128
        auto nParticles = 0;
        auto clusterEnergy = new Double_t[TA2GoAT_MAX_PARTICLE];
        treeTracks->Branch("nTracks", &nParticles, "nTracks/I");
        treeTracks->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nTracks]/D");

        // use some values to fill circularly
        const std::vector<double> possible_values{0, 23, 432, 12, 32, 6, 0.5}; // number is prime to increase randomness
        auto it_curr_value = possible_values.begin();

        for(auto fill=0;fill<nFills;fill++) {
            nParticles = fill;
            for(auto entry=0;entry<fill;entry++) {
                if(it_curr_value == possible_values.end())
                    it_curr_value = possible_values.begin();
                filled_data.push_back(*it_curr_value);
                clusterEnergy[entry] = *it_curr_value;
                ++it_curr_value;
            }
            treeTracks->Fill();
        }
    }

    REQUIRE(filled_data.size() == (nFills-1)*nFills/2);


}
