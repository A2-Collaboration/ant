#include "catch.hpp"
#include "catch_config.h"

#include "base/Interpolator.h"

#include "base/detail/interp2d/interp2d.h" // for INDEX_2D

#include <iostream>

using namespace std;
using namespace ant;

void dotest_symmetric(Interpolator2D::Type type);
void dotest_weird();

TEST_CASE("Interpolator2D: Bicubic", "[base]") {
    dotest_symmetric(Interpolator2D::Type::Bicubic);
}

TEST_CASE("Interpolator2D: Bilinear", "[base]") {
    dotest_symmetric(Interpolator2D::Type::Bilinear);
}


TEST_CASE("Interpolator2D: Weird stuff", "[base]") {
    dotest_weird();
}

void dotest_symmetric(Interpolator2D::Type type) {
    const vector<double> x{0.0, 1.0, 2.0, 3.0};
    const vector<double> y{0.0, 1.0, 2.0, 3.0};
    const vector<double> z{1.0, 1.1, 1.2, 1.3,
                           1.1, 1.2, 1.3, 1.4,
                           1.2, 1.3, 1.4, 1.5,
                           1.3, 1.4, 1.5, 1.6};
    Interpolator2D inter(x,y,z, type);

    // check with the same points
    for(size_t i=0;i<x.size();i++) {
        for(size_t j=0;j<y.size();j++) {
            CHECK(inter.GetPoint(x[i],y[j])
                  == Approx(z[INDEX_2D(i,j,x.size(),y.size())]));
        }
    }

    // check some extra points
    const vector<double> xval{0.0, 0.5, 1.0, 1.5, 2.5, 3.0};
    const vector<double> yval{0.0, 0.5, 1.0, 1.5, 2.5, 3.0};
    const vector<double> zval{1.0, 1.1, 1.2, 1.3, 1.5, 1.6};

    for(size_t i=0;i<xval.size();i++) {
        CHECK(inter.GetPoint(xval[i],yval[i]) == Approx(zval[i]));
    }
}

void dotest_weird() {
    REQUIRE_THROWS_AS(Interpolator2D inter({1,2,3,4},{1,2,3,4},{1,2,3}), Interpolator2D::Exception);
}



