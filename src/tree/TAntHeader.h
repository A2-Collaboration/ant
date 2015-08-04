#pragma once

#include "TNamed.h"

#include "tree/TDataRecord.h"

#include <string>

namespace ant {

class TAntHeader: public TNamed {
public:
    ant::TID FirstID;
    ant::TID LastID;
    std::string SetupName;
    std::string CmdLine;
    std::string WorkingDir;
    std::string GitInfo;

    TAntHeader(const std::string& title="");

    virtual ~TAntHeader();


    ClassDef(TAntHeader, ANT_UNPACKER_ROOT_VERSION)
};

}
