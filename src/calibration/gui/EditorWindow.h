#pragma once

#include "calibration/gui/Indicator_traits.h"
#include "calibration/Editor.h"

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

    class EditorCanvas : public TCanvas, public update_notify_traits
    {
        // update_notify_traits interface
    public:
        virtual void UpdateMe() override;
    };
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

    void CreateToolbar(TGVerticalFrame* frame);
    void UpdateLayout();

    MyComboBox* calibSelector;
    ant::calibration::Editor editor;

    bool flag;

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
