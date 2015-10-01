#ifndef ANT_SMARTTTREE
#define ANT_SMARTTTREE

#include <string>


class TTree;

namespace ant {

class SmartTree {
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



    virtual bool AddCut(const std::string& cut) =0;
    virtual bool RemoveCut(const std::string& cut) =0;



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

    static SmartTree* Create(TTree* tree);
};

}


#endif
