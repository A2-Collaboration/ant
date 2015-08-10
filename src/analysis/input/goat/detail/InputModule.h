#pragma once

#include <string>

class TTree;

namespace ant {
namespace analysis {
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
}
