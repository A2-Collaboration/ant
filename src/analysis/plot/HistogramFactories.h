#pragma once

#include "Histogram.h"
#include "SmartHist.h"
#include "tree/TParticle.h"

#include <string>

class TDirectory;
class TH1D;
class TH2D;
class TH3D;
class TTree;

namespace ant {
namespace analysis {

class SmartHistFactory {
private:

    TDirectory* dir = nullptr;

    TDirectory* goto_dir();
    void restore_dir(TDirectory* dir);

    HistogramFactory base_factory;

    std::string title_prefix;

    std::string MakeTitle(const std::string& title);

    /**
     * @brief create a new TDirectory. If a TDirectory already exists, append a number (1, 2, 3, ...) at the end
     * @param name Name of the directory to create
     * @param rootdir where to create
     * @return new created directory
     */
    static TDirectory* mkDirNumbered(const std::string& name, TDirectory* rootdir);


public:
    SmartHistFactory(const std::string& directory_name, TDirectory* root=nullptr, const std::string& title_prefix_ = "");
    SmartHistFactory(const std::string& directory_name, const SmartHistFactory &parent, const std::string& title_prefix_ = "");

    void SetRootDir(TDirectory* root_dir=nullptr);
    void SetTitlePrefix(const std::string& title_prefix_);
    void SetDirDescription(const std::string& desc);


    template<class Hist = TH1D>
    Hist* makeTH1D(
            const std::string& title,
            const std::string& xlabel,
            const std::string& ylabel,
            const BinSettings& bins,
            const std::string& name="")
    {

        TDirectory* old = goto_dir();
        Hist* r = base_factory.Make1D<Hist>(MakeTitle(title),xlabel,ylabel,bins,name);
        restore_dir(old);
        return r;
    }

    TH2D* makeTH2D(
            const std::string& title,
            const std::string& xlabel,
            const std::string& ylabel,
            const BinSettings& xbins,
            const BinSettings& ybins,
            const std::string& name="");

    TH3D* makeTH3D(const std::string& title,
            const std::string& xlabel,
            const std::string& ylabel,
            const std::string& zlabel,
            const BinSettings& xbins,
            const BinSettings& ybins,
            const BinSettings& zbins,
            const std::string& name="");

    TTree* makeTTree(const std::string& name);

    template<class T, typename... Args>
    T* make(Args&&... args) {
        TDirectory* old = goto_dir();
        auto t = new T(std::forward<Args>(args)...);
        restore_dir(old);
        return t;
    }

    // Predef smart hists

    SmartHist1<const TParticlePtr&> InvariantMass(
            const std::string& title,
            const std::string& xlabel="M [MeV]",
            const std::string& ylabel="",
            BinSettings bins=BinSettings(1000),
            const std::string& name="");

    SmartHist1<const TParticlePtr&> ThetaAngle(
            const std::string& title,
            const std::string& xlabel="#theta [#circ]",
            const std::string& ylabel="",
            BinSettings bins=BinSettings(180),
            const std::string& name="");

    SmartHist1<const TParticlePtr &> KinEnergyPlot(
            const std::string& title,
            const std::string& xlabel="#E_{k} [MeV]",
            const std::string& ylabel="",
            BinSettings bins=BinSettings(1000),
            const std::string& name="");


    template<typename T, typename FunctionType>
    SmartHist1<T> makeHist(FunctionType func,
            const std::string& title,
            const std::string& xlabel,
            const std::string& ylabel,
            const BinSettings& bins,
            const std::string& name="") {
        TDirectory* old = goto_dir();
            SmartHist1<T> r = SmartHist1<T>::makeHist(func,title,xlabel,ylabel,bins,name,base_factory);
        restore_dir(old);
        return std::move(r);
    }

    template<typename T>
    SmartHist1<T> makeHist(
            const std::string& title,
            const std::string& xlabel,
            const std::string& ylabel,
            const BinSettings& bins,
            const std::string& name="") {
        TDirectory* old = goto_dir();
            SmartHist1<T> r = SmartHist1<T>::makeHist(title,xlabel,ylabel,bins,name,base_factory);
        restore_dir(old);
        return std::move(r);
    }

    TH1D* copyTH1D(TH1D* hist, const std::string &newname);

    template<typename T>
    SmartHist1<T> Copy(const SmartHist1<T>& c, const std::string& newname) {
        TH1D* hist = copyTH1D(c.histogram, newname);
        return SmartHist1<T>(*hist, c.fillfunction->Copy());
    }

//    void Add(SmartHist1Base& hist, const string &as_name);
//    void Add(SmartHist1Base& hist);
};

}
}
