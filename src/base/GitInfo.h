#pragma once

#include <string>

namespace ant {

class GitInfo {

public:

    std::string GetUser() const;
    std::string GetDescription() const;
    bool IsDirty() const;


    GitInfo(const std::string& repofolder_ = "");
private:
    const std::string repofolder;
    std::string exec_git(const std::string& args) const;
};


}
