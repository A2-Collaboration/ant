#pragma once

#include "Rtypes.h"
#include "TNamed.h"

#include <string>
#include <vector>

#ifndef __CINT__
#include "analysis/plot/root_draw.h"
#endif

class THStack;
class TH1;

namespace ant {

/**
 * @brief The hstack class
 *
 * Wrapper for ROOT's THStack.
 * Overloads the << operator to add histograms to it, and of course
 * can be drawn the same way with the canvas class.
 * The axis labels are always taken from the last histogram added.
 * @see ant::canvas
 *
 */

class hstack : public TNamed
{

// for CINT, this class looks empty (except TNamed inheritance)
#ifndef __CINT__

protected:

    struct Hist_t {

        Hist_t(TH1* ptr, const std::string& option) : Ptr(ptr), Option(option) {}
        TH1* Ptr = nullptr;
        std::string Option;

        Hist_t() {}

        template<typename Archive>
        void load(Archive archive) {
            std::string path;
            archive(path, Option);
            Ptr = GetPtr(path);
        }
        template<typename Archive>
        void save(Archive archive) const {
            archive(GetPath(Ptr), Option);
        }

    protected:
        static TH1* GetPtr(const std::string& path);
        static std::string GetPath(const TH1* ptr);
    };


    std::vector<Hist_t> hists;

    std::string current_option;

    std::string xlabel;
    std::string ylabel;
    std::string title;

    bool UseIntelliLegend;
    bool IgnoreEmptyHist;

    void checkHists();

public:
    hstack(const std::string& name, const std::string& title="",
           bool useIntelliLegend = false,
           bool ignoreEmptyHist = false);


    hstack(hstack&&) = default;
    hstack& operator= (hstack&&) = default;

    hstack& operator<< (TH1* hist);
    hstack& operator<< (const drawoption& c);

    virtual void Draw(const char* option) override;
    virtual void Browse(TBrowser* b) override;

    template<typename Archive>
    void serialize(Archive archive) {
        archive(static_cast<TNamed&>(*this),
                hists, xlabel, ylabel, UseIntelliLegend, IgnoreEmptyHist);
        checkHists(); // see if any nullptr were loaded
    }

#endif

public:

    hstack();
    virtual ~hstack();
    ClassDef(hstack, 0)

private:
    // prevent ROOTcint from creating copy-constructors
    hstack(const hstack&);
    hstack& operator= (const hstack&);

};

}
