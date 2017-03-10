#pragma once

#include <string>
#include <map>

#include "mc/database/Query.h"


struct McAction {
    unsigned nEvents;
    std::string   outfile;
    double   Emin;
    double   Emax;
    virtual void Run() const =0;
    virtual ~McAction() = default;
};


using allowedTargets_t = std::map<std::string,ant::mc::data::Query::Selection>;
static const allowedTargets_t allowedTargets({
            {"proton",    ant::mc::data::Query::Selection::gpBeamTarget},
            {"neutron",   ant::mc::data::Query::Selection::gnBeamTarget}
        });

auto getAllowedTargetNames = [] ()
{
    std::vector<std::string> tn(allowedTargets.size());
    std::transform(allowedTargets.begin(),allowedTargets.end(),tn.begin(),
                   [](const std::pair<std::string,ant::mc::data::Query::Selection>& sel)
    {return sel.first;});
    return tn;
};
