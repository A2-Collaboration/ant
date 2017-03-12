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
     * @brief buildCmdLine little helper to reconstruct cmd line string from main(...) args
     * @param argc
     * @param argv
     * @return properly build string from argc and argv
     */
    static std::string buildCmdLine(int argc, const char* const* argv);

    /**
     * @brief list files in a directory
     * @param folder Path to directory to list
     * @param extension File extension to look for. all files if empty.
     * @return a list of all filenames with path prepended
     */
    static std::list<std::string> lsFiles(const std::string& folder=".",
                                          const std::string& extension="",
                                          bool ignoreDotDirs=false,
                                          bool doNotPrependFolder=false);

    /**
     * @brief testopen a given file
     * @param filename to open as test
     * @param errmsg non-empty if unsuccessful open
     * @return true if file can be opened for reading
     */
    static bool testopen(const std::string& filename, std::string& errmsg);

    static bool testopen(const std::string& filename);

    /**
     * @brief execute given cmd and return stdout (stderr is suppressed)
     * @param cmd the command to execute
     * @return stdout of cmd
     */
    static std::string exec(const std::string& cmd);

    /**
     * @brief getCwd determines current working directory (cwd)
     * @return properly constructed std::string containing cwd
     */
    static std::string getCwd();

    /**
     * @brief absolutePath returns the absolute path of relative path w.r.t. given cwd
     * @param path relative path
     * @param cwd if empty, the current working dir is used
     * @return absolute path
     */
    static std::string absolutePath(const std::string& path, const std::string& cwd = getCwd());

    /**
     * @brief Check if filename is a dead symlink
     * @param filename
     * @return
     */
    static bool isDeadLink(const std::string& filename);

    /**
     * @brief Check if a path/file exists. For symlinks: only check the link itself, independent from destination.
     * @param path
     * @return
     */
    static bool path_exists(const std::string& path);


};
}
}
