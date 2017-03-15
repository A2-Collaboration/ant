#pragma once
#include "analysis/plot/HistogramFactory.h"
#include "base/OptionsList.h"
#include "base/WrapTFile.h"
#include "base/std_ext/memory.h"
#include <functional>
#include <memory>
#include <vector>
#include <stdexcept>

namespace ant {
namespace analysis {


class Plotter {
private:
    std::string name_;

protected:
    HistogramFactory HistFac;

public:
    Plotter(const std::string& name, const WrapTFileInput& input, OptionsPtr opts);

    std::string GetName() const { return name_; }

    virtual long long GetNumEntries() const =0;
    virtual bool ProcessEntry(const long long entry) =0;
    virtual void Finish();
    virtual void ShowResult();

    virtual ~Plotter();

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

};

using plotter_creator = std::function< std::unique_ptr<Plotter>(const std::string& name, const WrapTFileInput& input, OptionsPtr otps) >;

class PlotterRegistry
{
    friend class PlotterRegistration;

private:
    using plotter_creators_t = std::map<std::string, plotter_creator>;
    plotter_creators_t plotter_creators;
    static PlotterRegistry& get_instance();

    void RegisterPlotter(plotter_creator c, const std::string& name) {
        plotter_creators[name] = c;
    }
public:

    static std::unique_ptr<Plotter> Create(const std::string& name, const WrapTFileInput &input, OptionsPtr opts = std::make_shared<OptionsList>());

    static std::vector<std::string> GetList();

    struct Exception : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

};

class PlotterRegistration
{
public:
    PlotterRegistration(plotter_creator c, const std::string& name);
};

template<class T>
std::unique_ptr<Plotter> plotter_factory(const std::string& name, const WrapTFileInput& input, OptionsPtr opts)
{
    return std_ext::make_unique<T>(name, input, opts);
}

#define AUTO_REGISTER_PLOTTER(plotter) \
    ant::analysis::PlotterRegistration _plotter_registration_ ## plotter(ant::analysis::plotter_factory<plotter>, #plotter);

}


}
