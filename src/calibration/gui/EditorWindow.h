#pragma once

#include "calibration/gui/Indicator_traits.h"
#include "calibration/Editor.h"

#include "TGFrame.h"
#include "KeySymbols.h"
#include "TCanvas.h"
#include "TGComboBox.h"
#include "TRootEmbeddedCanvas.h"

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

    class EditorCanvas : public TRootEmbeddedCanvas,
                           public update_notify_traits
    {
    protected:
        TCanvas* theCanvas;
    public:
        EditorCanvas(const TGWindow *p = 0) :
            TRootEmbeddedCanvas(0, p, 400, 400) // only important place to set some width/height
        {
            auto frame = (TGCompositeFrame*)fCanvasContainer;
            frame->RemoveInput(kKeyPressMask | kKeyReleaseMask);
            theCanvas = new TCanvas("Editor",10,10, GetCanvasWindowId());
            AdoptCanvas(theCanvas);
        }
        virtual void UpdateMe() override;

        virtual void cd(int subPad = 0){
            theCanvas->cd(subPad);
        }
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

    EditorCanvas* ecanvas;
    TGHorizontalFrame* frame_canvas = nullptr;

    ant::calibration::Editor editor;

    std::string currentCalID;

    TH2D* calHist;
    void drawCalibration();


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
