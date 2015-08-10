#include "root-addons.h"

#include "TCutG.h"
#include "TGraph.h"

using namespace ant::analysis::utils;

std::shared_ptr<TCutG> root::makeTCutG(const std::string& name, const std::initializer_list<std::pair<double, double> >& p) {

    auto cut = std::make_shared<TCutG>(name.c_str(), p.size());

    int i=0;
    for(auto& point : p) {
        cut->SetPoint(i++,point.first, point.second);
    }

    return move(cut);

}

std::shared_ptr<TGraph> root::makeTGraph(const std::initializer_list<std::pair<double, double> >& p) {

    auto cut = std::make_shared<TGraph>(p.size());

    int i=0;
    for(auto& point : p) {
        cut->SetPoint(i++,point.first, point.second);
    }

    return move(cut);

}
