#ifndef ANT_SMARTTTREE
#define ANT_SMARTTTREE

#include <string>

#include "TNamed.h"
#include "Rtypes.h"
#include "TCanvas.h"

class TTree;

namespace ant {

class SmartTreeCanvas : public TCanvas {
public:

    SmartTreeCanvas(const std::string& name, const std::string& title);
    virtual ~SmartTreeCanvas();

    virtual void SetFrozen(bool f=true) =0;  //*TOGGLE* *GETTER=GetFrozen
    virtual bool GetFrozen() const =0;

    virtual void SetCut(const char* cut) =0; //*MENU* *GETTER=GetCut
    virtual const char* GetCut() const =0;

    ClassDef(SmartTreeCanvas,1)
};

//class SmartTreeCanvas1D : public SmartTreeCanvas {
//public:
//    virtual void SetBinsX(UInt_t xbins) =0; //*MENU*
//};

class SmartTree: public TNamed {
protected:
    SmartTree(const std::string& name): TNamed(name.c_str(), "") {}

public:
    virtual ~SmartTree();

    virtual void Draw(const std::string& xepression, const int xbins=100) =0;
    virtual void Draw(const std::string& xexpression, const std::string& yexpression, const int xbins=100, const int ybins=100) =0;

    /**
     * @brief Set a range cut for an expression
     * @param expression
     * @param min
     * @param max
     */
    virtual void SetRange(const std::string& expression, const double min, const double max) =0;

    /**
     * @brief Remove range cut for an expression
     * @param expression The expression to remove to cut for
     */
    virtual void RemoveRange(const std::string& expression) =0;



    virtual bool Cut(const std::string& cut) =0;
    virtual bool RemoveCut(const std::string& cut) =0;
    virtual void UndoCut() =0;
    virtual void Limit(Long64_t n_entries=-1) =0;



    /**
     * @brief Print active cuts to stdout
     */
    virtual void PrintCuts() const =0;

    /**
     * @brief Remove all cuts.
     *        redraws all histograms of auto update is enabled, see SetAutoUpdate()
     */
    virtual void RemoveAllCuts() =0;

    /**
     * @brief Redraw all histograms and apply cuts
     */
    virtual void Update() =0;
    virtual bool GetAutoUpdate() const =0;
    virtual void SetAutoUpdate(bool update=true) =0;

    virtual void CloseAll() =0;

    static SmartTree* Create(TTree* tree);

    ClassDef(SmartTree, 1)
};

}


#endif
