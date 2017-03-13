#include "catch.hpp"

#include "base/WrapTTree.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/iterators.h"

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

    {
        // writing ROOTArray
        struct MyTree3 : WrapTTree {
            ADD_BRANCH_T(ROOTArray<double>, MyROOTArray)
        };
        MyTree3 t;
        TTree t_;
        REQUIRE_THROWS_AS(t.CreateBranches(addressof(t_)),WrapTTree::Exception);
    }
}

void dotest_pod() {
    constexpr auto nFills = 100; // not higher than max, see below

    tmpfile_t tmpfile;
    std::vector<std::vector<double>> filled_data;
    {
        WrapTFileOutput outputfile(tmpfile.filename, WrapTFileOutput::mode_t::recreate, true);
        // the following code is copied from
        // https://github.com/A2-Collaboration-dev/acqu/blob/master/acqu_user/root/src/TA2GoAT.h
        // https://github.com/A2-Collaboration-dev/acqu/blob/master/acqu_user/root/src/TA2GoAT.cc
        auto treeTracks = new TTree("treeTracks","treeTracks");
        constexpr auto TA2GoAT_MAX_PARTICLE = 128; // was #define TA2GoAT_MAX_PARTICLE 128
        auto nParticles = 0;
        auto clusterEnergy = new Double_t[TA2GoAT_MAX_PARTICLE];
        auto fbeam = new Float_t[5];
        treeTracks->Branch("nTracks", &nParticles, "nTracks/I");
        treeTracks->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nTracks]/D");
        // inspired by a2geant tree...
        treeTracks->Branch("beam", fbeam, "fbeam[5]/F");

        // use some values to fill circularly
        const std::vector<double> possible_values{3, 23, 432, 12, 32, 6, 0.5}; // number is prime to increase randomness
        auto it_value = std_ext::getCircularIterator(possible_values.begin(), possible_values.end());

        const std::vector<int> possible_sizes{5, 23, 12, 1, 128}; // number is another prime
        auto it_size = std_ext::getCircularIterator(possible_sizes.begin(), possible_sizes.end());

        for(auto fill=0;fill<nFills;fill++) {
            nParticles = *it_size;

            filled_data.emplace_back();

            fbeam[0] = *it_value;
            fbeam[1] = *it_value;
            fbeam[2] = *it_value;
            fbeam[3] = *it_value;
            fbeam[4] = *it_value;

            for(auto entry=0;entry<*it_size;entry++) {
                filled_data.back().push_back(*it_value);
                clusterEnergy[entry] = *it_value;
                ++it_value;
            }
            ++it_size;
            treeTracks->Fill();
        }

        delete[] clusterEnergy;
        // tree is owned by file, will close at this end-of-scope
    }

    REQUIRE_FALSE(filled_data.empty());

    // try reading the tree with WrapTTree::ROOTArray
    struct tree_t : WrapTTree {
        ADD_BRANCH_T(ROOTArray<Double_t>, clusterEnergy);
        ADD_BRANCH_T(ROOTArray<Float_t>,  beam);
    };

    tree_t t;

    WrapTFileInput inputfile(tmpfile.filename);
    REQUIRE(inputfile.GetObject("treeTracks", t.Tree));
    REQUIRE(t.Tree->GetEntries() == nFills);
    REQUIRE_NOTHROW(t.LinkBranches());

    for(auto entry=0; entry<t.Tree->GetEntries();entry++) {
        t.Tree->GetEntry(entry);
        REQUIRE(t.clusterEnergy().size()>0);
        REQUIRE(t.clusterEnergy() == filled_data[entry]);
        REQUIRE(t.beam().size()==5);
        REQUIRE(t.beam().at(0) == filled_data[entry][0]);
        REQUIRE(t.beam().at(1) == filled_data[entry][0]);
        REQUIRE(t.beam().at(2) == filled_data[entry][0]);
        REQUIRE(t.beam().at(3) == filled_data[entry][0]);
        REQUIRE(t.beam().at(4) == filled_data[entry][0]);
    }
}
