#include "unpacker/Unpacker.h"
#include "expconfig/ExpConfig.h"

#include "tree/THeaderInfo.h"
#include "base/Logger.h"
#include <iostream>

using namespace std;
using namespace ant;

int main(int argc, char** argv)
{
    SetupLogger();

    if(argc == 2) {

        auto unpacker = Unpacker::Get(argv[1]);
        LOG(INFO) << "Found setup: " <<  ExpConfig::Setup::GetLastFound()->GetName();

        while(auto item = unpacker->NextItem()) {
            auto headerInfo = dynamic_cast<THeaderInfo*>(item.get());
            if(headerInfo) {
                cout << *headerInfo << endl;
                return EXIT_SUCCESS;
            }
        }
        LOG(ERROR) << "Error: Did not find header info in file " << argv[1];
        return EXIT_FAILURE;
    }

    cout << "Try to figure out what file that is..." << endl;
    cout << "Usage: Ant-info <filename>" << endl;

    return EXIT_FAILURE;
}