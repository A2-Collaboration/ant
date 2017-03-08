#pragma once

#include <vector>
#include <string>
#include <cstdint>

namespace ant {

/**
 * @brief A tmpfolder class which deletes the whole directory on destroy
 */
struct tmpfolder_t {
    std::string foldername;
    tmpfolder_t();
    ~tmpfolder_t();

    // make it movable only, as the tmpfolder deletes its directory in dtor
    tmpfolder_t(tmpfolder_t&&) = default;
    tmpfolder_t(const tmpfolder_t&) = delete;
    tmpfolder_t& operator=(tmpfolder_t&&) = default;
    tmpfolder_t& operator=(const tmpfolder_t&) = delete;
};

/**
 * @brief A tmpfile class which cleans up after itself
 *
 * Mostly used by tests, but maybe helpful elsewhere
 */
struct tmpfile_t {
  std::string filename;
  std::vector<std::uint8_t> testdata;
  tmpfile_t();
  tmpfile_t(const tmpfolder_t& folder, const std::string& extension);
  void write_testdata();
  ~tmpfile_t();
  static std::size_t tmpfiles;

  // make it movable only, as the tmpfile deletes its file in dtor
  tmpfile_t(tmpfile_t&&) = default;
  tmpfile_t(const tmpfile_t&) = delete;
  tmpfile_t& operator=(tmpfile_t&&) = default;
  tmpfile_t& operator=(const tmpfile_t&) = delete;
};

}
