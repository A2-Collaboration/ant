#pragma once

#include <string>
#include "analysis/plot/Histogram.h"

class TTree;
class TH1;
class TH2;
class TH3;
class TCut;



TH1* Draw(TTree* tree, const std::string& formula, const TCut& cut, const int bins, const double min, const double max);
TH2* Draw(TTree* tree, const std::string& formula, const TCut& cut, const int xbins, const double xmin, const double xmax, const int ybins, const double ymin, const double ymax);

TH3* Draw(TTree* tree, const std::string& formula, const TCut& cut, const ant::analysis::BinSettings& xbins, const ant::analysis::BinSettings& ybins, const ant::analysis::BinSettings& zbins);
