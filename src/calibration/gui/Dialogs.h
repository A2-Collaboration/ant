#pragma once

#include <string>
#include <list>

#include "TGListBox.h"
#include "TGLabel.h"
#include "TExec.h"

namespace ant
{
namespace calibration
{
namespace gui{

class Query;

class DialogHandler_traits
{
friend class Query;
protected:
    virtual void captureReturnValue(Query* dialog, std::string& returnValue)=0;
};

class Query: public TExec
{
protected:
    DialogHandler_traits* owner;
    std::string returnValue;
public:
    Query(DialogHandler_traits* theOwner):
        owner(theOwner),
        returnValue("testRetVal")
    {}
//    void SendReturnValue() {owner->captureReturnValue(this,returnValue); }
    virtual void Exec(const char* arg) override
    {
        owner->captureReturnValue(this,returnValue);//     owner->SendReturnValue();
    }

};

class ListQuery: public Query
{
protected:
    TGMainFrame*       mainFrame;

    TGVerticalFrame*   vFrame;
    TGLabel*           label;
    TGListBox*         listBox;
    TGTextButton*      btnCancel;
    TGTextButton*      btnOk;

public:
    ListQuery(DialogHandler_traits* theOwner, const std::string& title, const std::string& text,const std::list<std::string>& items = {});

    void SetItemList(const std::list<std::string>& items);
    void ClearItemList();
    void AddItem(const std::string& item);
};

}
}
}
