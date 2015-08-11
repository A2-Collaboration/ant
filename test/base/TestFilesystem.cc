#include "catch.hpp"
#include "catch_config.h"

#include "base/tmpfile_t.h"
#include "base/filesystem.h"

#include <string>
#include <list>
#include <iostream>

using namespace std;
using namespace ant;

void dotest();
void makefiles();



TEST_CASE("Write WrapTfile", "[base]") {

    dotest();
}

void dotest() {
    tmpfile_t a(".tst");
    a.write_testdata();

    list<string> files;

    REQUIRE_NOTHROW(files = filesystem::lsFiles(".",".tst"));
    REQUIRE(files.size()==1);
    REQUIRE(files.front()==a.filename);

}

