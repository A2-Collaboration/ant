#pragma once

#include "Rtypes.h"
#include "THStack.h"

#include <string>
#include <vector>

#ifndef __CINT__
#include "analysis/plot/RootDraw.h"
#include "analysis/plot/HistStyle.h"
#include "base/interval.h"
#include "base/std_ext/math.h"
#include <memory>
#include <map>
#else
// disable the override keyword
#define override
#endif

class TH1;
class TLegend;
class TPaveText;

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

struct hstack_Menu;
struct hstack : THStack
{

// for CINT, this class looks empty (except TNamed inheritance and some methods)
#ifndef __CINT__

    using ModOption_t = analysis::plot::histstyle::ModOption_t;

    struct options_t {
        bool UseIntelliLegend = true ;
        bool IgnoreEmptyHist = true;
        bool ShowEntriesInLegend = true;
        bool UseIntelliTitle = true;
        bool FixLegendPosition = false;
        bool UseYAxisEntriesPerBin = false;
        interval<interval<double>> LegendPosition = {
            {0.7,  0.63},
            {0.88, 0.88}
        };
        interval<double> YAxisRange = {std_ext::NaN, std_ext::NaN};
        double MCScale = 1.0;
        struct hist_t {
            bool Hidden = false;
            double Scale = 1.0;
            std::string AddTo; // default empty
        };
        std::map<std::string, hist_t> PerHist;
        hist_t tryGetHist(const std::string& titlekey) const {
            return PerHist.find(titlekey) == PerHist.end() ? hist_t{} : PerHist.at(titlekey);
        }
        // ugly way of remembering the applied scaling
        struct scale_t {
            double Val = 1.0;
            operator double() { return Val; }
            scale_t& operator=(double v) { Val = v; return *this; }
        };

        std::map<TH1*, scale_t> AppliedScale;
    };

    hstack(const std::string& name, const std::string& title="", bool simple_ = false);

    bool IsCompatible(const hstack& other) const;

    bool Add(const hstack& other, const double scaling = 1);
    bool Divide(const hstack& other);
    bool Multiply(const hstack& other);

    hstack(hstack&&) = default;
    hstack& operator= (hstack&&) = default;

    hstack& operator<< (TH1* hist);
    hstack& operator<< (const drawoption& c);
    hstack& operator<< (const ModOption_t& option);


    friend std::ostream& operator<<( std::ostream& s, const hstack& o);

    template<typename Archive>
    void serialize(Archive& archive) {
        SetTitle(origtitle.c_str());
        archive(static_cast<TNamed&>(*this), hists, simple, cutnames);
        origtitle = GetTitle();
    }

    friend struct hstack_Menu;

    void UpdateMCScaling();

    static options_t GlobalOptions;

protected:

    using THStack::Add; // hide Add

    struct hist_t {

        hist_t(TH1* ptr, const ModOption_t& option) :
            Path(GetPath(ptr)),
            Ptr(ptr),
            Option(option)
        {}
        hist_t() {}

        std::string Path;
        TH1* Ptr = nullptr;
        ModOption_t Option;

        bool isDataHist() const;
        std::string getTitleKey() const;

        template<typename Archive>
        void load(Archive& archive) {
            archive(Path, Option);
            Ptr = GetPtr(Path);
        }
        template<typename Archive>
        void save(Archive& archive) const {
            archive(Path, Option);
        }

        static TH1* GetPtr(const std::string& path);
        static std::string GetPath(const TH1* ptr);
    };

    using hists_t = std::vector<hist_t>;
    hists_t hists;

    struct wraphist_t {
        wraphist_t(const hist_t& h) : Hist(h) {}
        // simple const ref does not work as it can't be swapped...
        std::reference_wrapper<const hist_t> Hist;
        double Entries = 0;
    };

    ModOption_t current_option;

    std::list<std::string> cutnames;
    std::string origtitle;
    bool simple;

    std::unique_ptr<TLegend>   intellilegend;
    std::unique_ptr<TPaveText> intellititle;

    void checkHists();
    void buildIntelliTitle();
    static void updateIntelliLegend(TLegend& legend, std::list<wraphist_t> wraphists);

#endif // __CINT__

public:

    virtual void Print(const char* option) const override;
    virtual void Print() const; // *MENU*
    virtual void Paint(const char* chopt) override;

    virtual void OpenStackMenu(); // *MENU*
    virtual void SetGlobalYAxisRange(double low, double high); // *MENU*

    virtual void FixLegendPosition(bool flag); // *TOGGLE* *GETTER=GetFixLegendPosition
    virtual bool GetFixLegendPosition() const;

    virtual void UseIntelliLegend(bool flag); // *TOGGLE* *GETTER=GetUseIntelliLegend
    virtual bool GetUseIntelliLegend() const;

    virtual void UseIntelliTitle(bool flag); // *TOGGLE* *GETTER=GetUseIntelliTitle
    virtual bool GetUseIntelliTitle() const;

    virtual void IgnoreEmptyHist(bool flag); // *TOGGLE* *GETTER=GetIgnoreEmptyHist
    virtual bool GetIgnoreEmptyHist() const;

    virtual void ShowEntriesInLegend(bool flag); // *TOGGLE* *GETTER=GetShowEntriesInLegend
    virtual bool GetShowEntriesInLegend() const;

    virtual void UseYAxisEntriesPerBin(bool flag); // *TOGGLE* *GETTER=GetUseYAxisEntriesPerBin
    virtual bool GetUseYAxisEntriesPerBin() const;

    // to be used with Ant-hadd
    virtual Long64_t Merge(TCollection* li, TFileMergeInfo *info) override;

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
