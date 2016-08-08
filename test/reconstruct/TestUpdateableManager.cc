#include "catch.hpp"

#include "reconstruct/UpdateableManager.h"
#include "reconstruct/Reconstruct_traits.h"

using namespace std;
using namespace ant;
using namespace ant::reconstruct;


void dotest1();
void dotest2();
void dotest3();
void dotest4();
void dotest5();
void dotest6();

TEST_CASE("UpdateableManager: Simple combinations", "[reconstruct]") {
    dotest1();
}

TEST_CASE("UpdateableManager: Update exactly at start", "[reconstruct]") {
    dotest2();
}

TEST_CASE("UpdateableManager: Start in between", "[reconstruct]") {
    dotest3();
}

TEST_CASE("UpdateableManager: Multiple updates", "[reconstruct]") {
    dotest4();
}

TEST_CASE("UpdateableManager: Unallowed things", "[reconstruct]") {
    dotest5();
}

TEST_CASE("UpdateableManager: Multiple updateables", "[reconstruct]") {
    dotest6();
}

// implement some testable Updateable item
struct UpdateableItem :  Updateable_traits {

    const list<TID> ChangePoints; // the available change points
    vector<TID> UpdatePoints;     // tracks the called updates


    UpdateableItem(const list<TID>& changePoints) :
        ChangePoints(changePoints)
    {}

    virtual std::list<Loader_t> GetLoaders() override
    {
        auto loader = [this] (const TID& currPoint, TID& nextChangePoint) {
            UpdatePoints.push_back(currPoint);
            for(auto& tid : ChangePoints) {
                if(tid > currPoint) {
                    nextChangePoint = tid;
                    break;
                }
            }
        };
        return {loader};
    }

};

// provide some points for testing
const vector<TID> p = {
    {0x00},
    {0x01},
    {0x10},
    {0x11},
    {0x12},
    {0x13}
};

void dotest1()
{
    // do some mocking with many points...
    for(const auto& startPoint : p) {
        for(const auto& updatePoint : p) {

            // item has no changepoints at all
            auto item = make_shared<UpdateableItem>(list<TID>{});

            UpdateableManager manager(startPoint, {item});
            manager.UpdateParameters(updatePoint);

            // check that the item was updated once
            REQUIRE(item->UpdatePoints.size() == 1);
        }
    }
}

void dotest2()
{
    auto item = make_shared<UpdateableItem>(list<TID>{p[1], p[3], p[5]});

    UpdateableManager manager(p[3], {item});
    manager.UpdateParameters(p[3]);

    REQUIRE(item->UpdatePoints.size() == 1);
    REQUIRE(p[3] == item->UpdatePoints.at(0));
}

void dotest3() {
    auto item = make_shared<UpdateableItem>(list<TID>{p[1], p[3], p[5]});

    UpdateableManager manager(p[2], {item});
    manager.UpdateParameters(p[4]);

    REQUIRE(item->UpdatePoints.size() == 2);
    REQUIRE(p[2] == item->UpdatePoints.at(0));
    REQUIRE(p[3] == item->UpdatePoints.at(1));
}

void dotest4() {
    auto item = make_shared<UpdateableItem>(list<TID>{p[1], p[3], p[5]});

    UpdateableManager manager(p[2], {item});

    manager.UpdateParameters(p[4]);
    manager.UpdateParameters(p[5]);

    REQUIRE(item->UpdatePoints.size() == 3);
    REQUIRE(p[2] == item->UpdatePoints.at(0));
    REQUIRE(p[3] == item->UpdatePoints.at(1));
    REQUIRE(p[5] == item->UpdatePoints.at(2));
}

void dotest5() {
    auto item = make_shared<UpdateableItem>(list<TID>{p[1], p[3], p[5]});

    UpdateableManager manager(p[2], {item});
    manager.UpdateParameters(p[4]);
    manager.UpdateParameters(p[4]);

    REQUIRE(item->UpdatePoints.size() == 2);
    REQUIRE(p[2] == item->UpdatePoints.at(0));
    REQUIRE(p[3] == item->UpdatePoints.at(1));
}

void dotest6() {
    auto item1 = make_shared<UpdateableItem>(list<TID>{p.begin(), p.end()});
    auto item2 = make_shared<UpdateableItem>(list<TID>{p[0], p[1], p[3]});

    UpdateableManager manager(p[0], {item1, item2});
    manager.UpdateParameters(p.back());

    REQUIRE(item1->UpdatePoints == p);

    REQUIRE(item2->UpdatePoints.size() == 3);

    REQUIRE(p[0] == item2->UpdatePoints.at(0));
    REQUIRE(p[1] == item2->UpdatePoints.at(1));
    REQUIRE(p[3] == item2->UpdatePoints.at(2));

}

