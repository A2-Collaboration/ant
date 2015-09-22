#pragma once

#include "calibration/Editor.h"
#include "base/interval.h"
#include "Dialogs.h"

#include <string>
#include <map>
#include <set>
#include <functional>
#include "Rtypes.h"
#include "analysis/plot/root_draw.h"
#include "calibration/gui/Indicator_traits.h"
#include "calibration/Editor.h"

#include "TCanvas.h"
#include "TPaveText.h"

#include "TRootEmbeddedCanvas.h"

class TH2D;


namespace ant {
namespace calibration {
namespace gui {

class EditorCanvas;

class EmbeddedEditorCanvas : public TRootEmbeddedCanvas,
        public update_notify_traits
{
private:
    EditorCanvas* theCanvas;

public:
    EmbeddedEditorCanvas(const std::shared_ptr<ant::calibration::Editor>& editor, const std::string& calID, const TGWindow *p = 0);

    virtual void SelectInvalid();
    virtual void SetCalID(const std::string& calID);
    virtual std::list<std::uint32_t> GetSelected();
    virtual void clearSelections();

    virtual void UpdateMe() override;
};


class EditorCanvas:
        public TCanvas,
        public update_notify_traits
{
private:


    TH2D*                         calHist;
    std::shared_ptr<Editor>       ed;

    std::string             currentCalID;

    std::set<std::uint32_t> indexMemory;
    interval<std::uint32_t> indexInterVal;
    bool                    intervalStartSet;

    void markInterval(Int_t y);
    void markLine(Int_t y);
    void updateCalHist();
    void HandleKeypress(const char key);

    void unFillLine(uint32_t lineNumber);
    void fillLine(uint32_t lineNumber);

public:
    EditorCanvas(const std::shared_ptr<Editor>& editor, const std::string& calID, int winID);
    void SetCalID(const std::string& calID);

    std::list<std::uint32_t> CreateSelectionList();

    virtual void UpdateMe() override;
    void ResetCalibration();
    void MarkInvalid();
    void HandleInput(EEventType button, Int_t x, Int_t y) override;
};

}
}
}
