#include "catch.hpp"
#include "catch_config.h"

#include "base/WrapTTree.h"
#include "base/tmpfile_t.h"
#include "base/WrapTFile.h"
#include "base/std_ext/memory.h"
#include "base/std_ext/iterators.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TChain.h"

using namespace std;
using namespace ant;

void dotest_basics();
void dotest_copy();
void dotest_nasty();
void dotest_pod();
void dotest_chain();
void dotest_opt_branches();
void dotest_stdarray();
void dotest_templating();


TEST_CASE("WrapTTree: Basics", "[base]") {
    dotest_basics();
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

TEST_CASE("WrapTTree: Read from TChain", "[base]") {
    dotest_chain();
}

TEST_CASE("WrapTTree: Optional branches", "[base]") {
    dotest_opt_branches();
}

TEST_CASE("WrapTTree: std::array", "[base]") {
    dotest_stdarray();
}

TEST_CASE("WrapTTree: templating", "[base]") {
    dotest_templating();
}


struct MyTree : WrapTTree {
    ADD_BRANCH_T(bool,           Flag1)        // simple type
    ADD_BRANCH_T(unsigned,       N1)           // simple type
    ADD_BRANCH_T(short,          N2)           // simple type
    ADD_BRANCH_T(bool,           Flag2, true)  // simple type
    ADD_BRANCH_T(vector<double>, SomeArray, 3) // use ctor of std::vector
    ADD_BRANCH_T(TLorentzVector, LV)           // complex ROOT type
};

void dotest_basics() {

    tmpfile_t tmpfile;

    // do some simple IO
    {
        WrapTFileOutput outputfile(tmpfile.filename);

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
        REQUIRE(t2.CopyFrom(t1));
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
        WrapTFileOutput outputfile(tmpfile.filename, true);
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
        ADD_BRANCH_T(ROOTArray<Double_t>, clusterEnergy)
        ADD_BRANCH_T(ROOTArray<Float_t>,  beam)
    };

    WrapTFileInput inputfile(tmpfile.filename);
    tree_t t;

    REQUIRE(inputfile.GetObject("treeTracks", t.Tree));
    REQUIRE(t.Tree->GetEntries() == nFills);
    REQUIRE_NOTHROW(t.LinkBranches());

    for(auto entry=0; entry<t.Tree->GetEntries();entry++) {
        t.Tree->GetEntry(entry);
        REQUIRE(t.clusterEnergy().size()>0);
        REQUIRE(t.clusterEnergy() == filled_data[entry]);
        REQUIRE(t.beam().size()==5);
        REQUIRE(t.beam[0] == filled_data[entry][0]);
        REQUIRE(t.beam[1] == filled_data[entry][0]);
        REQUIRE(t.beam[2] == filled_data[entry][0]);
        REQUIRE(t.beam[3] == filled_data[entry][0]);
        REQUIRE(t.beam[4] == filled_data[entry][0]);
    }
}


void dotest_chain() {

    auto make_treefile = [] (const std::string& filename) {
        WrapTFileOutput outputfile(filename, true);
        auto treeTracks = new TTree("tree","tree");
        auto nParticles = 0;
        auto clusterEnergy = new Double_t[128];
        auto beam = new Float_t[5];
        treeTracks->Branch("nTracks", &nParticles, "nTracks/I");
        treeTracks->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nTracks]/D");
        treeTracks->Branch("beam", beam, "beam[5]/F");

        for(auto entry=0;entry<10;entry++) {
            nParticles = 2;
            for(auto i=0;i<2;i++) {
                clusterEnergy[i] = 10*entry+i;
            }
            beam[0] = entry;
            beam[1] = entry;
            beam[2] = entry;
            beam[3] = entry;
            beam[4] = entry;
            treeTracks->Fill();
        }

        delete[] clusterEnergy;
    };

    // first, create a TChain of two TTree in file
    tmpfile_t tmpfile1;
    make_treefile(tmpfile1.filename);
    tmpfile_t tmpfile2;
    make_treefile(tmpfile2.filename);

    // try reading the tree with WrapTTree::ROOTArray
    struct tree_t : WrapTTree {
        ADD_BRANCH_T(ROOTArray<Double_t>, clusterEnergy)
        ADD_BRANCH_T(ROOTArray<Float_t>, beam)
    };

    // chain two files
    {
        auto chain = std_ext::make_unique<TChain>("tree");
        REQUIRE(chain->AddFile(tmpfile1.filename.c_str()) == 1);
        REQUIRE(chain->AddFile(tmpfile2.filename.c_str()) == 1);
        REQUIRE(chain->GetEntries() == 20);

        tree_t t;
        REQUIRE_NOTHROW(t.LinkBranches(chain.get()));

        for(int entry=0;entry<chain->GetEntries();entry++) {
            INFO("entry=" << entry);
            t.Tree->GetEntry(entry);
            const auto wrapped_entry = entry>9 ? entry-10 : entry;
            REQUIRE(t.clusterEnergy().size() == 2);
            REQUIRE(t.clusterEnergy[0] == 10*wrapped_entry+0);
            REQUIRE(t.clusterEnergy[1] == 10*wrapped_entry+1);
            const vector<float> expected(5, wrapped_entry);
            REQUIRE(t.beam() == expected);
        }
    }

    // chain just one file
    {
        auto chain = std_ext::make_unique<TChain>("tree");
        REQUIRE(chain->AddFile(tmpfile1.filename.c_str()) == 1);
        REQUIRE(chain->GetEntries() == 10);

        tree_t t;
        REQUIRE_NOTHROW(t.LinkBranches(chain.get()));

        for(int entry=0;entry<chain->GetEntries();entry++) {
            INFO("entry=" << entry);
            t.Tree->GetEntry(entry);
            REQUIRE(t.clusterEnergy().size() == 2);
            REQUIRE(t.clusterEnergy[0] == 10*entry+0);
            REQUIRE(t.clusterEnergy[1] == 10*entry+1);
            const vector<float> expected(5, entry);
            REQUIRE(t.beam() == expected);
        }
    }
}

void dotest_opt_branches() {
    struct MyOptTree_t : WrapTTree {
        ADD_BRANCH_T(double, SomeNumber)
        ADD_BRANCH_T(std::vector<int>, SomeIndices)
        ADD_BRANCH_T(bool, SomeFlag)
        ADD_BRANCH_OPT_T(bool, SomeOptionalFlag)
        ADD_BRANCH_OPT_T(std::vector<double>, SomeOptionalNumbers)
    };

    auto create_opt_tree = [] (const string& filename, bool skipOptional) {
        WrapTFileOutput outputfile(filename);
        MyOptTree_t t;
        REQUIRE_FALSE(t.SomeOptionalFlag.IsPresent);
        REQUIRE_FALSE(t.SomeOptionalNumbers.IsPresent);
        t.CreateBranches(outputfile.CreateInside<TTree>("test","test"), skipOptional);
        t.SomeFlag = true;
        t.SomeIndices().push_back(1);
        t.SomeIndices().push_back(2);
        t.SomeIndices().push_back(3);
        t.SomeIndices().push_back(5);
        t.SomeOptionalFlag = true;
        t.SomeOptionalNumbers().resize(4);
        t.SomeOptionalNumbers[0] = 5;
        t.SomeOptionalNumbers[1] = 4;
        t.SomeOptionalNumbers[2] = 3;
        t.SomeOptionalNumbers[3] = 2;

        t.Tree->Fill();
    };

    // do not skip optional branches on creation
    {
        tmpfile_t tmpfile;
        create_opt_tree(tmpfile.filename, false);

        WrapTFileInput inputfile(tmpfile.filename);
        MyOptTree_t t;
        REQUIRE_FALSE(t.SomeOptionalFlag.IsPresent);
        REQUIRE_FALSE(t.SomeOptionalFlag);

        REQUIRE(inputfile.GetObject("test", t.Tree));
        t.LinkBranches();
        REQUIRE(t.SomeOptionalFlag.IsPresent);
        REQUIRE(t.SomeOptionalNumbers.IsPresent);

        t.Tree->GetEntry();
        bool myFlag = t.SomeOptionalFlag;
        std::vector<double> optionalNumbers(t.SomeOptionalNumbers);
        REQUIRE(myFlag);
        const vector<double> expected{5,4,3,2};
        REQUIRE(optionalNumbers == expected);
    }

    // skip optional branches on creation
    {
        tmpfile_t tmpfile;
        create_opt_tree(tmpfile.filename, true);

        WrapTFileInput inputfile(tmpfile.filename);
        {
            MyOptTree_t t;
            REQUIRE_FALSE(t.SomeOptionalFlag.IsPresent);
            REQUIRE_FALSE(t.SomeOptionalFlag);

            REQUIRE(inputfile.GetObject("test", t.Tree));
            REQUIRE_NOTHROW(t.LinkBranches());
            REQUIRE_FALSE(t.SomeOptionalFlag.IsPresent);
            REQUIRE_FALSE(t.SomeOptionalNumbers.IsPresent);

            t.Tree->GetEntry();
            REQUIRE_FALSE(t.SomeOptionalFlag);
            REQUIRE(t.SomeOptionalNumbers().empty());
        }
        {
            MyOptTree_t t;
            REQUIRE(inputfile.GetObject("test", t.Tree));
            REQUIRE_THROWS_AS(t.LinkBranches(true), WrapTTree::Exception);
        }
    }

}

void dotest_stdarray() {

    // make two-dimensional branch (inspired by 'dircos' of a2geant h12 tree)
    tmpfile_t tmpfile;
    {
        WrapTFileOutput outputfile(tmpfile.filename, true);
        auto tree = new TTree("tree","tree");
        auto npart = 0;
        Float_t dircos[100][3];
        tree->Branch("npart", &npart, "npart/I");
        tree->Branch("dircos", dircos, "dircos[npart][3]/F");
        tree->Branch("otherI", dircos, "other1[npart][5]/I");
        tree->Branch("otherD", dircos, "other1[npart][5]/D");

        dircos[npart][0] = 1;
        dircos[npart][1] = 2;
        dircos[npart++][2] = 3;
        dircos[npart][0] = 4;
        dircos[npart][1] = 5;
        dircos[npart++][2] = 6;
        tree->Fill();

        npart = 0;
        dircos[npart][0] = 7;
        dircos[npart][1] = 8;
        dircos[npart++][2] = 9;
        tree->Fill();

    }

    WrapTFileInput inputfile(tmpfile.filename);

    struct MyStdArrayTree : WrapTTree {
        ADD_BRANCH_T(ROOTArray_Float<3>, dircos)
    };
    MyStdArrayTree t;
    REQUIRE(inputfile.GetObject("tree", t.Tree));
    REQUIRE(t.Tree->GetEntries()==2);

    // linking dircos is not supported (yet)!
    REQUIRE_NOTHROW(t.LinkBranches());

    {
        t.Tree->GetEntry(0);
        REQUIRE(t.dircos().size() == 2);
        const std::vector<std::array<float,3>> expected{
            {{1.0f,2.0f,3.0f}}, {{4.0f,5.0f,6.0f}}
        };
        REQUIRE(t.dircos() == expected);
    }
    {
        t.Tree->GetEntry(1);
        REQUIRE(t.dircos().size() == 1);
        const std::vector<std::array<float,3>> expected{
            {{7.0f,8.0f,9.0f}}
        };
        REQUIRE(t.dircos() == expected);
    }

}

template<typename MyTree = WrapTTree>
class MyClass {

    struct Tree_t : MyTree {
        ADD_BRANCH_T(double, SomeBranch)
    };
};

void dotest_templating() {
    MyClass<> test;
}
