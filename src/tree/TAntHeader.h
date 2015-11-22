#pragma once

#include "TID.h"

#include "TNamed.h"

#include <string>

namespace ant {

/**
 * @brief The TAntHeader class holds information about the histogram file
 */
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

    /**
     * @brief Print AntHeader inforamtion
     * @note This is an override of the ROOT method. All options are ignored.
     * @see Print()
     */
    virtual void Print(Option_t*) const;

    /**
     * @brief Print AntHeader information
     * @note Does the same as the other print method. This is jut for convenience (no arguments).
     */
    virtual void Print() const; //*MENU*


    ClassDef(TAntHeader, ANT_UNPACKER_ROOT_VERSION)
};

}
