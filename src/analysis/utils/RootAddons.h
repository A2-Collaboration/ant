#pragma once

#include <memory>
#include <string>

class TCutG;
class TGraph;

namespace ant {
namespace analysis {
namespace utils {

/**
 * @brief The root struct
 */
struct root {
/**
 * @brief Create a TCutG from a list of points
 * @param name The name to give the new object.
 * @param p list of points of doubles (x,y). Example {{1,1},{2,3},{4,2}}
 * @return a new TGutG as shared pointer
 * @bug Don't give the same name, including "" to more than one TCutG. Might lead to crashes during ROOT cleanup. This is a ROOT bug.
 */
static std::shared_ptr<TCutG> makeTCutG(const std::string& name, const std::initializer_list<std::pair<double,double>>& p);

/**
 * @brief Create a TGraph from a list of points
 * @param name The name to give the new object.
 * @param p list of points of doubles (x,y). Example {{1,1},{2,3},{4,2}}
 * @return a new TGraph as shared pointer
  */
static std::shared_ptr<TGraph> makeTGraph(const std::initializer_list<std::pair<double,double>>& p);

};

}
}
}
