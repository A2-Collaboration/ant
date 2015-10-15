#pragma once

#include <string>
#include <list>
#include <vector>
#include <map>
#include <memory>
#include <functional>

#include "Rtypes.h"

class THStack;
class TH1D;
class TCanvas;
class TObject;
class TVirtualPad;

namespace ant {

class canvas;

class root_drawable_traits {
    friend class ant::canvas;
public:
    virtual void Draw(const std::string& option) const =0;
    virtual ~root_drawable_traits() = default;
};

class modifier {
public:
  virtual ~modifier() = default;
};

class endcanvas: public modifier {
public:
    endcanvas() {}
};

static const endcanvas endc;

struct rowend : modifier {
    rowend() {}
};

static const rowend endr;

class drawoption: public modifier {
protected:
    std::string option;
public:
    drawoption(const std::string& opt=""): option(opt) {}
    virtual ~drawoption() {}
    const std::string& Option() const { return option; }
};

enum class padoption_t {
    LogX, LogY, LogZ,
    Legend
};

class padoption: public modifier {

public:

    using map_options_t = std::map<padoption_t, std::function<void(TVirtualPad*)> >;
    static const map_options_t map_options;

protected:
    class base {
    protected:
        padoption_t option;
    public:
        base(const padoption_t& _option) : option(_option) {}

        const padoption_t& Option() const { return option; }
    };

public:

    class set : public base {
    public:
        set(const padoption_t& option) : base(option) {}
    };

    class unset : public base {
    public:
        unset(const padoption_t& option) : base(option) {}
    };

    virtual ~padoption() {}

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

    TCanvas* create(const std::string& title="");
    TCanvas* find();

    template<typename T>
    class drawable_container: public ant::root_drawable_traits {
        public:
            T obj;

            drawable_container(T o): obj(o) {}
            void Draw(const std::string &option) const {
                obj->Draw(option.c_str());
            }
        drawable_container(const drawable_container&) = delete;
        drawable_container& operator= (const drawable_container&) = delete;
    };

    using ObjectOption =
        std::tuple<
            std::unique_ptr<root_drawable_traits>,
            std::string,              // draw option
            std::list<padoption_t>,    // pad options
            bool // end of row
        >;

    std::list<ObjectOption> objs;

    std::string current_drawoption;
    std::list<padoption_t> current_padoptions;

    static bool isMarker(const std::unique_ptr<root_drawable_traits>& t);

    void DrawObjs(TCanvas* c, unsigned cols, unsigned rows);
public:

    canvas(const std::string& title="");
    virtual ~canvas();

    virtual void cd();

    virtual canvas& operator<< (root_drawable_traits& drawable);

    virtual canvas& operator<< (TObject* hist);

    virtual canvas& operator<< (const endcanvas& c);

    virtual canvas& operator<< (const rowend&);

    virtual canvas &operator<< (const drawoption& c);

    virtual canvas &operator<< (const padoption::set& c);

    virtual canvas &operator<< (const padoption::unset& c);

    virtual canvas &operator>> (const std::string& filename);

};

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
class hstack: public root_drawable_traits {
protected:
    THStack* stack;
    std::string current_option;
    std::string xlabel;
    std::string ylabel;

public:
    hstack(const std::string& name, const std::string &title="");
    virtual ~hstack();

    hstack(const hstack&) = delete;
    hstack& operator= (const hstack&) = delete;

    virtual hstack &operator<< (TH1D* hist);
    virtual hstack &operator<< (const drawoption& c);

    void Draw(const std::string &option) const override;
};

struct ColorPalette {
    static const std::vector<Color_t> Colors;
};

}
