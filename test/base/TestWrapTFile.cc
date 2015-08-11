#include "catch.hpp"
#include "catch_config.h"

#include "base/WrapTFile.h"
#include "base/std_ext.h"

#include <cstdio>
#include <memory>

#include "TCutG.h"


using namespace std;
using namespace ant;

void writeFile(const std::string& filename);
void readFile(const std::string& filename);

TEST_CASE("Write WrapTfile", "[base]") {
    writeFile("anttmp.root");
}

TEST_CASE("Read WrapTfile", "[base]") {
    readFile("anttmp.root");
    remove("anttmp.root");
    ///@todo Need a tmpfile_t that creates file names that end in .root
}

void writeFile(const std::string& filename) {

    TCutG* cut = new TCutG("test",3);
    cut->SetPoint(0,1,1);
    cut->SetPoint(1,2,1);
    cut->SetPoint(2,3,1);
    cut->SaveAs(filename.c_str()); // filename must end in .root
}

/**
 * @brief readFile test: shows that we can read objects from a root file and close the file, object still valid
 * @param filename File to open
 */
void readFile(const std::string& filename) {
    std::unique_ptr<WrapTFileInput> file = ant::std_ext::make_unique<ant::WrapTFileInput>();
    std::shared_ptr<TCutG> cut;

    REQUIRE_NOTHROW( file->OpenFile(filename));
    REQUIRE_NOTHROW( cut = file->GetSharedClone<TCutG>("test"));
    REQUIRE(cut != nullptr);
    REQUIRE_NOTHROW(file=nullptr);
    REQUIRE(cut->GetN() == 3);
}

