#pragma once

#include <string>

#include "TGListBox.h"
#include "TGLabel.h"

namespace ant
{
namespace calibration
{

class ListQueryWindow
{
protected:
    TGMainFrame*       mainFrame;

    TGVerticalFrame*   vFrame;
    TGLabel*           label;
    TGListBox*         listBox;

    TGHorizontalFrame* hFrame;
    TGTextButton*      btnCancel;
    TGTextButton*      btnOk;

public:
    ListQueryWindow(const std::string& title, const std::string& text,const std::list<std::string>& items)
    {
        mainFrame = new TGMainFrame(gClient->GetRoot(),10,10,kVerticalFrame);


        vFrame    = new TGVerticalFrame(mainFrame);
        mainFrame->AddFrame(vFrame,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

        label     = new TGLabel(vFrame,text.c_str());
        vFrame->AddFrame(label,new TGLayoutHints(kLHintsCenterX | kLHintsTop));

        listBox   = new TGListBox(vFrame);
        vFrame->AddFrame(listBox,new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX|kLHintsExpandY));

        int i = 0;
        for (const auto& item: items )
            listBox->AddEntry(item.c_str(),++i);


        hFrame    = new TGHorizontalFrame(mainFrame);
        mainFrame->AddFrame(vFrame,new TGLayoutHints(kLHintsExpandX | kLHintsBottom));

        btnCancel = new TGTextButton(hFrame,"&btnCancel");
        hFrame->AddFrame(btnCancel,new TGLayoutHints(kLHintsExpandX,5,5,3,4));

        btnOk     = new TGTextButton(hFrame,"&btnOk");
        hFrame->AddFrame(btnOk,new TGLayoutHints(kLHintsExpandX,5,5,3,4));

        mainFrame->SetCleanup(kDeepCleanup);



        mainFrame->MapSubwindows();
        mainFrame->Resize();

        mainFrame->SetWMSizeHints(mainFrame->GetDefaultWidth(),
                                  mainFrame->GetDefaultHeight(),
                                  1000,1000,0,0);
        mainFrame->SetWindowName(title.c_str());
        mainFrame->MapRaised();
    }


};

}
}
