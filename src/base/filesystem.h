#pragma once

#include <string>
#include <list>

namespace ant {

/**
 * @brief The filesystem struct contains static functions to work with files and folders.
 * @note This should be a namespace. But somehow this does not work... so struct it is (fix?)
 */
struct filesystem {

    /**
     * @brief list files in a directory
     * @param path Path to directory to list
     * @param extension File extension to look for. all files if empty.
     * @return a list of all filenames
     */
    static std::list<std::string> lsFiles(const std::string& path=".", const std::string& extension="");
};

}
