#pragma once

#include <string>
#include <list>
#include <vector>
#include <memory>
#include <functional>
#include "base/BinSettings.h"

#include "Rtypes.h"

class TCanvas;
class TObject;
class TVirtualPad;
class TTree;

namespace ant {

class canvas;

struct root_drawable_traits {
    virtual void Draw(const std::string& option) const =0;
    virtual ~root_drawable_traits() = default;
};

struct TTree_drawable : root_drawable_traits {

    // Basic constructor, use the same as drawing from TTree
    TTree_drawable(TTree* tree, const std::string& formula, const std::string& cut = "");

    // 1D constructor, using same form as Histogram factory (overloaded so that if no name is given, title is used)
    TTree_drawable(TTree* tree, const std::string& varx, const std::string& cut, const std::string &title, const std::string &xlabel, const std::string &ylabel, const BinSettings &binsx, const std::string &name);
    TTree_drawable(TTree* tree, const std::string& varx, const std::string& cut, const std::string &title, const std::string &xlabel, const std::string &ylabel, const BinSettings &binsx) :
        TTree_drawable(tree,varx,cut,title,xlabel,ylabel,binsx,title){}

    // 2D constructor, using same form as Histogram factory (overloaded so that if no name is given, title is used)
    TTree_drawable(TTree* tree, const std::string& varx, const std::string& vary, const std::string& cut, const std::string &title, const std::string &xlabel, const std::string &ylabel, const BinSettings &binsx, const BinSettings &binsy, const std::string &name);
    TTree_drawable(TTree* tree, const std::string& varx, const std::string& vary, const std::string& cut, const std::string &title, const std::string &xlabel, const std::string &ylabel, const BinSettings &binsx, const BinSettings &binsy) :
        TTree_drawable(tree,varx,vary,cut,title,xlabel,ylabel,binsx,binsy,title){}

    virtual void Draw(const std::string& option) const override;
    static unsigned nInstances;
protected:
    TTree* Tree;
    std::string Formula;
    std::string Cut;
    std::string Title;
    std::string Name;
    std::string Xlabel;
    std::string Ylabel;
    bool autolabels;
};

struct endcanvas {};

extern const endcanvas endc;

struct endrow {};

extern const endrow endr;

struct samepad_t {};

extern const samepad_t samepad;

struct drawoption {
    drawoption(const std::string& opt=""): option(opt) {}
    const std::string& Option() const { return option; }
protected:
    std::string option;
};

using padmodifier_t = std::function<void(TVirtualPad*)>;

struct padoption : padmodifier_t {

protected:
    struct permanent {
        permanent(const padmodifier_t& m) : Modifier(m) {}
        const padmodifier_t Modifier;
    };

    using padmodifier_t::padmodifier_t;

public:

    static const padoption Legend;
    static const padoption LogX;
    static const padoption LogY;
    static const padoption LogZ;
    static const padoption MakeSquare;

    struct SetFillColor : padmodifier_t {
        explicit SetFillColor(Color_t col);
    };

    struct enable : permanent {
        using permanent::permanent;
    };

    struct disable : permanent {
        using permanent::permanent;
    };
};

/**
 * @brief The canvas class
 *
 * Wrapper for TCanvas.
 * Drawable classes, either TObject or classes inhereting from root_drawable_traits,
 * can be added via the << operator. Just like with cout:
 * canvas c("test");
 * c << hist1 << hist2 << endc;
 *
 * The canvas is automatically subdivided to fit all added histograms.
 * To finally draw the canvas add the endc obejct to it.
 *
 * Draw options can be set using the drawoption modifier:
 *
 * c << drawoption("colz") << hist2d_a << hist2d_b << endc;
 */
class canvas {
protected:
    static unsigned int num;

    const std::string name;
    const std::string title;
    bool automode = true;
    bool addobject = false;

    TCanvas* FindTCanvas() const;

    struct DrawableItem {
        std::shared_ptr<root_drawable_traits> Drawable;
        std::string Option;
        DrawableItem(std::shared_ptr<root_drawable_traits> drawable,
                     const std::string& option) :
            Drawable(move(drawable)),
            Option(option)
        {}
    };


    struct pad_t  {
        std::list<DrawableItem> DrawableItems;     // objects to draw on pad (empty indicates end row)
        using PadOptions_t = std::list<padmodifier_t>;
        PadOptions_t PadOptions; // pad options
        pad_t(const PadOptions_t& options) : PadOptions(options) {}
        pad_t() {}
    };

    std::list<pad_t> pads;

    std::string current_drawoption;
    std::list<padmodifier_t> global_padoptions;
    std::list<padmodifier_t> onetime_padoptions;


    void DrawObjs(TCanvas* c, unsigned cols, unsigned rows);
    void AddDrawable(std::shared_ptr<root_drawable_traits> drawable);

public:

    canvas(const std::string& title_="");
    ~canvas();

    void cd();
    void clear();

    template<typename Drawable>
    // use SFINAE/enable_if to restrict this templated operator to types deriving from root_drawable_traits only
    // the return type of this operator<< is simply "canvas&" as for all others
    typename std::enable_if<std::is_base_of<root_drawable_traits, Drawable>::value, canvas>::type&
    operator<< (const Drawable& drawable) {
        // the type Drawable must be copyable, we use that here
        AddDrawable(std::make_shared<Drawable>(drawable));
        return *this;
    }

    canvas& operator<< (const std::shared_ptr<root_drawable_traits>& drawable) {
        AddDrawable(drawable);
        return *this;
    }

    canvas& operator<< (TObject* hist);

    canvas& operator<< (const endcanvas&);

    canvas& operator<< (const endrow&);

    canvas& operator<< (const samepad_t&);

    canvas& operator<< (const drawoption& c);

    canvas& operator<< (const padmodifier_t& c);

    canvas& operator<< (const padoption::enable& c);

    canvas& operator<< (const padoption::disable& c);

    canvas& operator>> (const std::string& filename);

};

struct ColorPalette {
    static const std::vector<Color_t> Colors;
};



}
