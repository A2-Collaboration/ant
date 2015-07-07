#ifndef INPUTMODULE_H
#define INPUTMODULE_H

#include "FileManager.h"

class TTree;

namespace ant {
namespace input {

class TreeRequestManager {
public:
    virtual ~TreeRequestManager () {}
    virtual TTree* GetTree(const std::string& name) =0;
};

class BaseInputModule {
public:
    virtual ~BaseInputModule() {}

    virtual bool SetupBranches(TreeRequestManager&& input_files) =0;
    virtual void GetEntry() =0;

};
}
}

#endif
