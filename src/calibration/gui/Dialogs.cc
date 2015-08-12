#include "Dialogs.h"


using namespace std;
using namespace ant;
using namespace ant::calibration::gui;



void ListQuery::AddItem(const string& item)
{
    entryStrings.emplace_back(item);
    listBox->AddEntry(item.c_str(), listBox->GetNumberOfEntries());
}

void ListQuery::ClearItemList()
{
    entryStrings.clear();
    listBox->RemoveAll();
}

void ListQuery::SetItemList(const list<string>& items)
{
    ClearItemList();
    for (const auto& item: items)
        AddItem(item);
}

void ListQuery::SetReturnValue()
{
    if (returnValue.size() == 0)
        returnValue.push_back(entryStrings.at(listBox->GetSelected()));
    else returnValue.at(0) =  entryStrings.at(listBox->GetSelected());
}

ListQuery::ListQuery(DialogHandler_traits* theOwner, const string& title, const string& text, const list<string>& items):
    Query(theOwner)
{

    mainFrame = new TGMainFrame(gClient->GetRoot(),10,10,kVerticalFrame);


    vFrame    = new TGVerticalFrame(mainFrame);
    mainFrame->AddFrame(vFrame,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    label     = new TGLabel(vFrame,text.c_str());
    vFrame->AddFrame(label,new TGLayoutHints(kLHintsCenterX | kLHintsTop));

    listBox   = new TGListBox(vFrame);
    vFrame->AddFrame(listBox,new TGLayoutHints(kLHintsTop|kLHintsLeft|kLHintsExpandX|kLHintsExpandY));

   // btnCancel = new TGTextButton(vFrame,"&Cancel");
  //  vFrame->AddFrame(btnCancel,new TGLayoutHints(kLHintsExpandX|kLHintsBottom|kLHintsLeft
 //                                                ,5,5,3,4));
//    btnCancel->Connect("Pressed()","TExec",this,"Exec(\"exit\")");

    btnOk = new TGTextButton(vFrame,"&Select");
    vFrame->AddFrame(btnOk,new TGLayoutHints(kLHintsExpandX|kLHintsBottom|kLHintsRight
                                             ,5,5,3,4));
    btnOk->Connect("Pressed()","TExec",this,"Exec(\"\")");

    for (const auto& item: items )
        AddItem(item);


    mainFrame->SetCleanup(kDeepCleanup);



    mainFrame->MapSubwindows();
    mainFrame->Resize();

    mainFrame->SetWMSizeHints(mainFrame->GetDefaultWidth(),
                              mainFrame->GetDefaultHeight(),
                              1000,1000,0,0);
    mainFrame->SetWindowName(title.c_str());
    mainFrame->MapRaised();
}

#include <iostream>

void Query::Exec(const char* signal)
{
    std::string the_signal(signal);
    cout << signal << endl;
    if (the_signal.find("exit") == 0)
    {
        this->Delete();
        return;
    }
    SetReturnValue();
    owner->captureReturnValue(this,returnValue);
}
