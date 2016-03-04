#pragma once

#include "Rtypes.h"
#include "TNamed.h"

#include <string>
#include <vector>

#ifndef __CINT__
#include "analysis/plot/root_draw.h"
#include "analysis/plot/HistStyle.h"
#include "base/printable.h"
#include "base/interval.h"
#else
// disable the override keyword
#define override
#endif

class THStack;
class TH1;

namespace ant {

/**
 * @brief The hstack class
 *
 * Intelligent wrapper for ROOT's THStack.
 *
 * Overloads the << operator to add histograms to it, and of course
 * can be drawn the same way with the canvas class.
 *
 * The axis labels are always taken from the last histogram added.
 * @see ant::canvas
 *
 * Extra features include intelligent legend building,
 * it can be written to ROOT files (provided the referenced TH1's are also present),
 * it can be merged (with Ant-hadd in order to load the class ant::hstack properly)
 *
 */

#ifndef __CINT__
struct hstack :  TNamed, printable_traits
#else
struct hstack : TNamed
#endif
{

// for CINT, this class looks empty (except TNamed inheritance and some methods)
#ifndef __CINT__

    using ModOption_t = analysis::plot::histstyle::ModOption_t;

    struct options_t {
        bool UseIntelliLegend = true ;
        bool IgnoreEmptyHist = true;
        bool ShowEntriesInLegend = true;
        bool UseIntelliTitle = true;
    };

    hstack(const std::string& name, const std::string& title="");

    bool IsCompatible(const hstack& other) const;

    hstack(hstack&&) = default;
    hstack& operator= (hstack&&) = default;

    hstack& operator<< (TH1* hist);
    hstack& operator<< (const drawoption& c);
    hstack& operator<< (const ModOption_t& option);


    virtual std::ostream& Print( std::ostream& s) const override;

    template<typename Archive>
    void serialize(Archive archive) {
        archive(static_cast<TNamed&>(*this),
                hists, xlabel, ylabel);
        checkHists();
    }

    static double Global_MC_Scaling;
    static std::map<TH1*, double> Scaled_Hists;
    void UpdateMCScaling();

    static interval<double> GlobalYAxisRange;
    static interval<interval<double>> GlobalLegendPosition;

    static options_t GlobalOptions;

protected:


    struct hist_t {

        hist_t(TH1* ptr, const ModOption_t& option) :
            Path(GetPath(ptr)),
            Ptr(ptr),
            Option(option)
        {}

        // clear the Ptr on copy
        // checkHists() needs to be called then
        hist_t(const hist_t& other) {
            Path = other.Path;
            Option = other.Option;
            Ptr = nullptr;
        }
        hist_t& operator= (const hist_t&) = delete;

        // move is ok
        hist_t& operator= (hist_t&&) = default;
        hist_t(hist_t&&) = default;

        std::string Path;
        TH1* Ptr = nullptr;
        ModOption_t Option;

        hist_t() {}

        bool isDataHist() const;

        template<typename Archive>
        void load(Archive archive) {
            archive(Path, Option);
            Ptr = GetPtr(Path);
        }
        template<typename Archive>
        void save(Archive archive) const {
            archive(Path, Option);
        }

        static TH1* GetPtr(const std::string& path);
        static std::string GetPath(const TH1* ptr);
    };

    using hists_t = std::vector<hist_t>;
    hists_t hists;

    ModOption_t current_option;

    std::string xlabel;
    std::string ylabel;
    std::string title;


    void checkHists();
    void buildIntelliLegend() const;
    void buildIntelliTitle() const;

#endif // __CINT__

public:

    virtual void Print(const char* option) const override;
    virtual void Print() const; // *MENU*
    virtual void Draw(const char* option) override;
    virtual void Browse(TBrowser* b) override;

    virtual void SetGlobalMCScaling(double scaling); // *MENU*
    virtual void SetGlobalYAxisRange(double low, double high); // *MENU*
    virtual void SetGlobalLegendPosition(double x1=0.5, double y1=0.67, double x2=0.88, double y2=0.88); // *MENU*

    virtual void UseIntelliLegend(bool flag); // *TOGGLE* *GETTER=GetUseIntelliLegend
    virtual bool GetUseIntelliLegend() const;

    virtual void UseIntelliTitle(bool flag); // *TOGGLE* *GETTER=GetUseIntelliTitle
    virtual bool GetUseIntelliTitle() const;

    virtual void IgnoreEmptyHist(bool flag); // *TOGGLE* *GETTER=GetIgnoreEmptyHist
    virtual bool GetIgnoreEmptyHist() const;

    virtual void ShowEntriesInLegend(bool flag); // *TOGGLE* *GETTER=GetShowEntriesInLegend
    virtual bool GetShowEntriesInLegend() const;

    // to be used with Ant-hadd
    Long64_t Merge(TCollection* li);

    hstack();
    virtual ~hstack();
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#endif
    ClassDef(hstack, 0)
#ifdef __clang__
#pragma clang diagnostic pop
#endif

private:
    // prevent ROOTcint from creating copy-constructors
    hstack(const hstack&);
    hstack& operator= (const hstack&);

};

}
