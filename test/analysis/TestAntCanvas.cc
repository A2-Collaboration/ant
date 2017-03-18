#include "catch.hpp"

#include "analysis/plot/RootDraw.h"

#include "TTree.h"
#include "TCanvas.h"
#include "TList.h"

#include <iostream>
#include <set>

using namespace std;
using namespace ant;

struct drawable_tester : root_drawable_traits {
    mutable int drawCalls = 0;

    virtual void Draw(const string&) const override {
        ++drawCalls;
    }
};

struct canvas_tester : canvas {
    using canvas::canvas;

    struct grid_t {
        size_t sizeX;
        size_t sizeY;
    };

    grid_t ReconstructGrid() {
        auto c = FindTCanvas();

        auto lnk = (TObjOptLink*)c->GetListOfPrimitives()->FirstLink();
        vector<TVirtualPad*> pads;
        while (lnk) {
            auto pad = dynamic_cast<TVirtualPad*>(lnk->GetObject());
            if(pad)
                pads.push_back(pad);
            lnk = (TObjOptLink*)lnk->Next();
        }

        // it is evil to do that (as doubles are compared), but it works here
        std::set<double> x;
        std::set<double> y;

        for(auto& pad : pads) {
            x.insert(pad->GetAbsXlowNDC());
            y.insert(pad->GetAbsYlowNDC());
        }

        return {x.size(), y.size()};
    }
};

TEST_CASE("AntCanvas: Simple AutoGrids", "[analysis]") {

    // 1 drawable
    {
        auto drawable = make_shared<drawable_tester>();
        canvas_tester c("test1");
        c << drawable
          << endc;
        auto grid = c.ReconstructGrid();
        REQUIRE(grid.sizeX == 1);
        REQUIRE(grid.sizeY == 1);
        REQUIRE(drawable->drawCalls == 1);
    }

    // 2 drawable
    {
        auto drawable = make_shared<drawable_tester>();
        canvas_tester c("test1");
        c << drawable
          << drawable
          << endc;
        auto grid = c.ReconstructGrid();
        REQUIRE(grid.sizeX == 2);
        REQUIRE(grid.sizeY == 1);
        REQUIRE(drawable->drawCalls == 2);
    }

    // 3 drawable
    {
        auto drawable = make_shared<drawable_tester>();
        canvas_tester c("test1");
        c << drawable
          << drawable
          << drawable
          << endc;
        auto grid = c.ReconstructGrid();
        REQUIRE(grid.sizeX == 2);
        REQUIRE(grid.sizeY == 2);
        REQUIRE(drawable->drawCalls == 3);
    }

    // 5 drawable
    {
        auto drawable = make_shared<drawable_tester>();
        canvas_tester c("test1");
        c << drawable
          << drawable
          << drawable
          << drawable
          << drawable
          << endc;
        auto grid = c.ReconstructGrid();
        REQUIRE(grid.sizeX == 3);
        REQUIRE(grid.sizeY == 2);
        REQUIRE(drawable->drawCalls == 5);
    }
}

TEST_CASE("AntCanvas: Manual ant::endr", "[analysis]") {

    // one endr
    {
        auto drawable = make_shared<drawable_tester>();
        canvas_tester c("test1");
        c << drawable
          << endr
          << drawable
          << endc;
        auto grid = c.ReconstructGrid();
        REQUIRE(grid.sizeX == 1);
        REQUIRE(grid.sizeY == 2);
        REQUIRE(drawable->drawCalls == 2);
    }

    // one endr, with one endr superfluous
    {
        auto drawable = make_shared<drawable_tester>();
        canvas_tester c("test1");
        c << drawable
          << endr
          << drawable
          << endr
          << endc;
        auto grid = c.ReconstructGrid();
        REQUIRE(grid.sizeX == 1);
        REQUIRE(grid.sizeY == 2);
        REQUIRE(drawable->drawCalls == 2);
    }

    // two endr, with one endr superfluous
    {
        auto drawable = make_shared<drawable_tester>();
        canvas_tester c("test1");
        c << drawable << drawable
          << endr
          << drawable << drawable
          << endr
          << endc;
        auto grid = c.ReconstructGrid();
        REQUIRE(grid.sizeX == 2);
        REQUIRE(grid.sizeY == 2);
        REQUIRE(drawable->drawCalls == 4);
    }

}