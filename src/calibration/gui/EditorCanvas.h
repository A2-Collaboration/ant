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
class TCalibrationData;
namespace calibration {
namespace gui {

class EditorCanvas;
class EditorWindow;

class EmbeddedEditorCanvas : public TRootEmbeddedCanvas,
        public update_notify_traits
{
private:
    EditorCanvas* theCanvas;

public:
    EmbeddedEditorCanvas(EditorWindow* EditorWindow,  const TGWindow *p = 0);

    virtual void EditSelection();
    virtual void SetToAverage();

    virtual void UpdateMe() override;
};


class EditorCanvas:
        public TCanvas,
        public update_notify_traits
{
private:

    TH1D*                         calDataHist;

    EditorWindow*       editorWindow;
    std::shared_ptr<ant::calibration::Editor> editor;

    void applyDataChanges();

    void HandleKeypress(const char key);

public:
    EditorCanvas(EditorWindow* EditorWindow, int winID);

    virtual void UpdateMe() override;

//    void ResetCalibration();
    void StartEditData();
    void SetToAverage();

    void HandleInput(EEventType button, Int_t x, Int_t y) override;
};

}
}
}
