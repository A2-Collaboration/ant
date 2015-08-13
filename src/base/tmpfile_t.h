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

struct tmpfolder_t {
    std::string foldername;
    tmpfolder_t();
    ~tmpfolder_t();
};

struct tmpfile_t {
  std::string filename;
  std::vector<std::uint8_t> testdata;
  tmpfile_t();
  tmpfile_t(const tmpfolder_t& folder, const std::string& extension);
  void write_testdata();
  ~tmpfile_t();
  static std::size_t tmpfiles;
};



}
