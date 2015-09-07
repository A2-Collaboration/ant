#pragma once

#include <list>
#include <string>

namespace ant {
namespace std_ext {

struct system {

    /**
     * @brief Check if running an interactive shell / if a TTY is attached to stdin
     * @return
     */
    static bool isInteractive();

    /**
     * @brief list files in a directory
     * @param folder Path to directory to list
     * @param extension File extension to look for. all files if empty.
     * @return a list of all filenames with path prepended
     */
    static std::list<std::string> lsFiles(const std::string& folder=".", const std::string& extension="");

};
}
}
