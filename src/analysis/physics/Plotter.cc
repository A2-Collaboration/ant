#include "Plotter.h"

using namespace ant;
using namespace ant::analysis;
using namespace std;

PlotterRegistry& PlotterRegistry::get_instance()
{
    static PlotterRegistry instance;
    return instance;
}

std::unique_ptr<Plotter> PlotterRegistry::Create(const string &name, const WrapTFileInput& input, OptionsPtr opts)
{
    auto& creators = PlotterRegistry::get_instance().plotter_creators;

    auto creator = creators.find(name);
    if(creator == creators.end())
        throw std::runtime_error("Plotter class " + name + " not found");

    // this may throw an exception
    auto plotter = creator->second(name, input, opts);

    return plotter;
}

std::vector<string> PlotterRegistry::GetList()
{
    std::vector<std::string> list;
    for(const auto& entry : get_instance().plotter_creators) {
        list.emplace_back(entry.first);
    }
    return list;
}


PlotterRegistration::PlotterRegistration(plotter_creator c, const string& name)
{
    PlotterRegistry::get_instance().RegisterPlotter(c,name);
}

Plotter::Plotter(const string &name, const WrapTFileInput&, OptionsPtr):
    name_(name),
    HistFac(name)
{}

void Plotter::Finish() {}

void Plotter::ShowResult() {}

ant::analysis::Plotter::~Plotter() {}
