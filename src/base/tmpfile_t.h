#pragma once

#include <vector>
#include <string>
#include <cstdint>


/**
 * @brief A tmpfile class which cleans up after itself
 *
 * Mostly used by tests, but maybe helpful elsewhere
 */

namespace ant {

/**
 * @brief A self-deleting tempfile.
 */
struct tmpfile_t {
    std::string filename;
    std::vector<std::uint8_t> testdata;

    /**
     * @brief Construtor. Generates a temp filename of the from anttmpfile.XXXXXX" + extension
     * @param extension optional file extension.
     */
    tmpfile_t(const std::string& extension="");
    void write_testdata();
    ~tmpfile_t();
};

}
