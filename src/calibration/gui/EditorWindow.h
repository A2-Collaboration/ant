#pragma once

#include "TGFrame.h"
#include "KeySymbols.h"

#include <functional>
#include <list>
#include <map>


namespace ant {
namespace calibration {
namespace gui {

class EditorWindow : public TGMainFrame
{
private:



    std::map<EKeySym, TGTextButton*> keys;


    void CreateToolbar(TGVerticalFrame* frame);
    void UpdateLayout();

    bool flag;

public:
    EditorWindow();
    virtual Bool_t HandleKey(Event_t *event) override;
    virtual ~EditorWindow();

    EditorWindow(const EditorWindow&) = delete;
    EditorWindow& operator=(const EditorWindow&) = delete;
};
}
}
}
