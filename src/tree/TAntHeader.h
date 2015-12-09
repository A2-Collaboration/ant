#pragma once

#include "TID.h"

#include "TNamed.h"

#include <string>

namespace ant {

/**
 * @brief The TAntHeader class holds information about the histogram file
 */
#ifndef __CINT__
struct TAntHeader : TNamed, printable_traits
#else
struct TAntHeader : TNamed
#endif
{
    ant::TID FirstID;
    ant::TID LastID;
    std::string SetupName;
    std::string CmdLine;
    std::string WorkingDir;
    std::string GitInfo;
    std::string GitInfoDatabase;

    TAntHeader(const std::string& title="");

    /**
     * @brief IsCompatible checks if other TAntHeader might be from same calibration iteration
     * @param other
     * @return
     */
    bool IsCompatible(const TAntHeader& other) const;


#ifndef __CINT__
    virtual std::ostream& Print( std::ostream& s) const override;
#endif

    // for convenience when used within ROOT shell
    virtual void Print(Option_t*) const;
    virtual void Print() const; //*MENU*

    virtual ~TAntHeader();
    ClassDef(TAntHeader, ANT_UNPACKER_ROOT_VERSION)
};

}
