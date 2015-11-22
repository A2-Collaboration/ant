#include "TAntHeader.h"

#include <iostream>

#include "TDirectory.h"

using namespace std;

ant::TAntHeader::TAntHeader(const std::string& title):
    TNamed("AntHeader", title.c_str()) {
    gDirectory->Add(this);
}

ant::TAntHeader::~TAntHeader()
{
}

void ant::TAntHeader::Print(Option_t*) const
{
    cout << "TAntHeader:\n"
         << "  FirstID:   " << FirstID << "\n"
         << "  LastID:    " << LastID  << "\n"
         << "  Setup:     " << SetupName << "\n"
         << "  CmdLine:   " << CmdLine << "\n"
         << "  WorkingDir:" << WorkingDir << "\n"
         << "  GitInfo:   " << GitInfo << endl;
}

void ant::TAntHeader::Print() const
{
    Print("");
}
