#include "catch.hpp"

#include "base/TH_ext.h"
#include "base/std_ext/string.h"

#include "TH3.h"
#include "TRandom.h"

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Test TH_ext GetSlice","[base]")
{
    dotest();
}

void dotest() {
    TH3D h("h", "random hist", 20,-10,10, 20,-10,10, 20,0,20);

    for (int i = 0; i < 1000; i++)
        h.Fill(gRandom->Gaus(-2,5),
               gRandom->Gaus(3,4),
               i%20);

    string info = "Test projection: ";
    for (const char* projection : {"yx", "xy", "xz"}) {  // test different projections (remember to reset range if necessary!)
        INFO(info + projection);
        if (std_ext::contains(string(projection), "xz"))  // reset z axis range to initial values
            h.GetZaxis()->SetRange(1, 20);
        for (int bin : {3, 7, 12, 15, 18}) {  // test slice on different bins per projection
            TH2D* slice = TH_ext::GetSlice(h, bin, projection);
            REQUIRE(slice != nullptr);
            // use SetRange on the correct axis
            if (std_ext::contains(string(projection), "xz"))
                h.GetYaxis()->SetRange(bin, bin);
            else
                h.GetZaxis()->SetRange(bin, bin);
            TH2D* proj = dynamic_cast<TH2D*>(h.Project3D(projection));

            bool match = true;
            for (int binx = 0; binx <= h.GetNbinsX()+1; ++binx)
                for (int biny = 0; biny <= h.GetNbinsY()+1; ++biny)
                    if (slice->GetBinContent(slice->GetBin(binx,biny))
                            != proj->GetBinContent(proj->GetBin(binx,biny)))
                        match = false;

            REQUIRE(match);
            REQUIRE(slice->GetRMS() == proj->GetRMS());
            REQUIRE(slice->GetMean(2) == proj->GetMean(2));
            REQUIRE(slice->GetEntries() == proj->GetEntries());
        }
    }
}
