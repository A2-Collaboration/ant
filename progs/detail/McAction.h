#pragma once

#include <string>

struct McAction {
    unsigned nEvents;
    std::string   outfile;
    double   Emin;
    double   Emax;
    virtual void Run() const =0;
    virtual ~McAction() = default;
};
