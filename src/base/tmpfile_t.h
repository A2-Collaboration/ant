#ifndef TMPFILE_T_H
#define TMPFILE_T_H

#include <vector>
#include <string>
#include <cstdint>


/**
 * @brief A tmpfile class which cleans up after itself
 *
 * Mostly used by tests, but maybe helpful elsewhere
 */

namespace ant {

struct tmpfile_t {
  std::string filename;
  std::vector<std::uint8_t> testdata;
  tmpfile_t();
  void write_testdata();
  ~tmpfile_t();
};

}

#endif // TMPFILE_T_H
