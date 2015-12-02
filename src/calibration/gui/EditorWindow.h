#pragma once

#include "calibration/Editor.h"
#include "calibration/gui/Indicator_traits.h"

#include "TGFrame.h"
#include "KeySymbols.h"
#include "TCanvas.h"
#include "TGComboBox.h"

#include <functional>
#include <vector>
#include <list>
#include <map>
#include <string>



namespace ant {
namespace calibration {
namespace gui {

class EmbeddedEditorCanvas;
class EditorCanvas;

class EditorWindow : public TGMainFrame,
        public update_notify_traits
{
private:

    std::shared_ptr<ant::calibration::Editor> editor;
    std::map<EKeySym, TGTextButton*> keys;

    TGTextButton* rootButton_avg;
    TGTextButton* rootButton_saveQuit;

    void createToolbar(TGVerticalFrame* frame);
    void updateLayout();


    EmbeddedEditorCanvas* ecanvas;
    TGHorizontalFrame* frame_canvas = nullptr;



    TH2D* calHist;

public:
    EditorWindow(const std::string& filename);

    virtual Bool_t HandleKey(Event_t *event) override;

    virtual ~EditorWindow();
    EditorWindow(const EditorWindow&) = delete;
    EditorWindow& operator=(const EditorWindow&) = delete;

    std::shared_ptr<ant::calibration::Editor> GetEditor();

    // update_notify_traits interface
public:
    void UpdateMe() override;
};


}
}
}
