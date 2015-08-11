#pragma once

#include <string>
#include <list>

#include "TGListBox.h"
#include "TGLabel.h"

namespace ant
{
namespace calibration
{
namespace gui{

class ListQuery
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
    ListQuery(const std::string& title, const std::string& text,const std::list<std::string>& items);
};

}
}
}
