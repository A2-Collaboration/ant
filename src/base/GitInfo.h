#pragma once

#include <string>

namespace ant {

class GitInfo {
public:
//    static GitInfo& I() {
//        // Returns the only instance
//        // Guaranteed to be lazy initialized
//        // Guaranteed that it will be destroyed correctly
//        static GitInfo instance;
//        return instance;
//    }

//    // forbid copies
//    GitInfo(const GitInfo&) = delete;
//    GitInfo& operator=(const GitInfo&) = delete;

    static std::string GetUser();
    static std::string GetDescription();

    // use static functions only
    GitInfo() = delete;
private:
    static std::string exec_git(const std::string& args);
};


}
