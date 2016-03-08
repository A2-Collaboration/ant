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
    virtual void captureReturnValue(Query* dialog, std::vector<std::string>& returnValue)=0;
protected:
    ~DialogHandler_traits() = default;
};

class Query: public TExec
{
protected:
    DialogHandler_traits* owner;
    std::vector<std::string> returnValue;
    virtual void SetReturnValue()=0;
public:
    Query(DialogHandler_traits* theOwner):
        owner(theOwner),
        returnValue()
    {}
//    void SendReturnValue() {owner->captureReturnValue(this,returnValue); }
    virtual void Exec(const char* signal) override;
};


class ListQuery: public Query
{
protected:
    TGMainFrame*       mainFrame;
    TGVerticalFrame*   vFrame;
    TGListBox*         listBox;
    TGLabel*           label;
//    TGTextButton*      btnCancel;
    TGTextButton*      btnOk;

    std::vector<std::string> entryStrings;
    void SetReturnValue() override;

public:
    ListQuery(DialogHandler_traits* theOwner, const std::string& title, const std::string& text,const std::list<std::string>& items = {});

    void SetItemList(const std::list<std::string>& items);
    void ClearItemList();
    void AddItem(const std::string& item);

    // Query interface
};

}
}
}
