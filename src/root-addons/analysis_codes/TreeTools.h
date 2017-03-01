#pragma once

#include <string>
#include "analysis/plot/HistogramFactory.h"

class TTree;
class TH1;
class TH2;
class TH3;
class TCut;



TH1* Draw(TTree* tree, const std::string& formula, const TCut& cut, const std::string& xtitle, const std::string& ytitle, const ant::BinSettings& xbins, const std::string& name);
TH2* Draw(TTree* tree, const std::string& formula, const TCut& cut, const std::string& xtitle, const std::string& ytitle, const ant::BinSettings& xbins, const ant::BinSettings& ybins, const std::string& name);

TH3* Draw(TTree* tree, const std::string& formula, const TCut& cut, const ant::BinSettings& xbins, const ant::BinSettings& ybins, const ant::BinSettings& zbins);
