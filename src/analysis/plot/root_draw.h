#pragma once

#include <string>
#include <list>
#include <vector>
#include <map>
#include <memory>
#include <functional>

#include "Rtypes.h"

class THStack;
class TH1;
class TH1D;
class TCanvas;
class TObject;
class TVirtualPad;

namespace ant {

class canvas;

struct root_drawable_traits {
    virtual void Draw(const std::string& option) const =0;
    virtual ~root_drawable_traits() = default;
};

struct canvas_modifier {};

struct endcanvas : canvas_modifier {};

extern const endcanvas endc;

struct endrow : canvas_modifier {};

extern const endrow endr;

struct samepad_t : canvas_modifier {};

extern const samepad_t samepad;

struct drawoption : canvas_modifier {
    drawoption(const std::string& opt=""): option(opt) {}
    const std::string& Option() const { return option; }
protected:
    std::string option;
};

using padmodifier_t = std::function<void(TVirtualPad*)>;

struct padoption : canvas_modifier, padmodifier_t {

protected:
    struct permanent {
        permanent(const padmodifier_t& m) : Modifier(std::addressof(m)) {}
        const padmodifier_t* Modifier;
    };

    using padmodifier_t::padmodifier_t;

public:

    static const padoption Legend;
    static const padoption LogX;
    static const padoption LogY;
    static const padoption LogZ;


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

    std::string name;
    bool automode = true;
    bool addobject = false;

    TCanvas* CreateTCanvas(const std::string& title="");
    TCanvas* FindTCanvas();

    template<typename T>
    struct drawable_container :  root_drawable_traits {
        T Object;

        drawable_container(T object) :
            Object(object)  {}

        void Draw(const std::string& option) const {
            Object->Draw(option.c_str());
        }

        drawable_container(const drawable_container&) = delete;
        drawable_container& operator= (const drawable_container&) = delete;
    };

    struct DrawableItem {
        std::unique_ptr<root_drawable_traits> Drawable;
        std::string Option;
        DrawableItem(std::unique_ptr<root_drawable_traits> drawable,
                     const std::string& option) :
            Drawable(move(drawable)),
            Option(option)
        {}
    };


    struct pad_t  {
        std::list<DrawableItem> DrawableItems;     // objects to draw on pad (empty indicates end row)
        using PadOptions_t = std::list<const padmodifier_t*>;
        PadOptions_t PadOptions; // pad options
        pad_t(const PadOptions_t& options) : PadOptions(options) {}
        pad_t() {}
    };

    std::list<pad_t> pads;

    std::string current_drawoption;
    std::list<const padmodifier_t*> global_padoptions;
    std::list<const padmodifier_t*> onetime_padoptions;


    void DrawObjs(TCanvas* c, unsigned cols, unsigned rows);
    void AddDrawable(std::unique_ptr<root_drawable_traits> drawable);

public:

    canvas(const std::string& title="");
    ~canvas();

    void cd();

    canvas& operator<< (root_drawable_traits& drawable);

    canvas& operator<< (TObject* hist);

    canvas& operator<< (const endcanvas& c);

    canvas& operator<< (const endrow&);

    canvas& operator<< (const samepad_t&);

    canvas& operator<< (const drawoption& c);

    canvas& operator<< (const padoption& c);

    canvas& operator<< (const padoption::enable& c);

    canvas& operator<< (const padoption::disable& c);

    canvas& operator>> (const std::string& filename);

};

struct ColorPalette {
    static const std::vector<Color_t> Colors;
};

}
