#pragma once

#include <string>

namespace ant {

class GitInfo {

public:

    std::string GetUser();
    std::string GetDescription();

    GitInfo(const std::string& repofolder_ = "");
private:
    const std::string repofolder;
    std::string exec_git(const std::string& args);
};


}
