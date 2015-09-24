#pragma once

#include "calibration/Editor.h"
#include "calibration/gui/EditorCanvas.h"
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

class EditorWindow : public TGMainFrame
{
private:

    class MyComboBox : public TGComboBox
    {
    public:
        MyComboBox(const TGWindow *p = 0, Int_t id = -1,
              UInt_t options = kHorizontalFrame | kSunkenFrame | kDoubleBorder,
              Pixel_t back = GetWhitePixel());
        void SetList(const std::list<std::string>& items);
        std::string GetSelectedText();
    };

    std::map<EKeySym, TGTextButton*> keys;

    TGTextButton* rootButton_delete;
    TGTextButton* rootButton_markInValid;
    TGTextButton* rootButton_StartEditor;

    MyComboBox* calibSelector;

    void createSelector(TGVerticalFrame* frame);
    void createToolbar(TGVerticalFrame* frame);
    void updateLayout();


    EmbeddedEditorCanvas* ecanvas;
    TGHorizontalFrame* frame_canvas = nullptr;

    std::shared_ptr<ant::calibration::Editor> editor;

    std::string dataFolder;
    std::string currentCalID;

    TH2D* calHist;

    void deleteSelections();

    void disableButtons();
public:
    EditorWindow(const std::string& folder);

    virtual Bool_t HandleKey(Event_t *event) override;

    virtual ~EditorWindow();
    EditorWindow(const EditorWindow&) = delete;
    EditorWindow& operator=(const EditorWindow&) = delete;
};
}
}
}
