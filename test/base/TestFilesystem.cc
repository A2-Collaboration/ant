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
    tmpfolder_t folder;
    tmpfile_t a(folder, ".tst");
    a.write_testdata();

    list<string> files;

    REQUIRE_NOTHROW(files = filesystem::lsFiles(folder.foldername,".tst"));
    REQUIRE(files.size()==1);
    REQUIRE(files.front()==string(a.filename));


    tmpfile_t b(folder, ".tst");
    b.write_testdata();

    REQUIRE_NOTHROW(files = filesystem::lsFiles(folder.foldername,".tst"));
    REQUIRE(files.size()==2);
}

