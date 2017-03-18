#pragma once

#include "calibration/Editor.h"

#include <functional>
#include "Rtypes.h"
#include "analysis/plot/RootDraw.h"
#include "calibration/gui/Indicator_traits.h"

#include "TCanvas.h"
#include "TRootEmbeddedCanvas.h"

class TH1D;


namespace ant {

struct TCalibrationData;

namespace calibration {

struct Editor;

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
    virtual void ResetData();
    virtual void ApplyChanges();

    virtual void UpdateMe() override;
};


class EditorCanvas:
        public TCanvas,
        public update_notify_traits
{
private:

    TH1D*          calDataHist;
    EditorWindow*  editorWindow;

    std::shared_ptr<ant::calibration::Editor> editor;


    void HandleKeypress(const char key);

public:
    EditorCanvas(EditorWindow* EditorWindow, int winID);

    virtual void UpdateMe() override;

    void ApplyDataChanges();
    void ResetCalibration();
    void StartEditData();
    void SetToAverage();

    void HandleInput(EEventType button, Int_t x, Int_t y) override;
};

}
}
}
