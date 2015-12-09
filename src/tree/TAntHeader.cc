#include "TAntHeader.h"

#include <iostream>

#include "TDirectory.h"


using namespace std;
using namespace ant;

TAntHeader::TAntHeader(const std::string& title):
    TNamed("AntHeader", title.c_str()) {
    gDirectory->Add(this);
}

TAntHeader::~TAntHeader()
{

}

bool TAntHeader::IsCompatible(const TAntHeader& other) const {
    return SetupName == other.SetupName
            && GitInfo == other.GitInfo
            && GitInfoDatabase == other.GitInfoDatabase;
}

ostream& TAntHeader::Print(ostream& s) const {
    s << "TAntHeader:\n"
      << "  FirstID:         " << FirstID << "\n"
      << "  LastID:          " << LastID  << "\n"
      << "  Setup:           " << SetupName << "\n"
      << "  CmdLine:         " << CmdLine << "\n"
      << "  WorkingDir:      " << WorkingDir << "\n"
      << "  GitInfo:         " << GitInfo << "\n"
      << "  GitInfoDatabase: " << GitInfoDatabase << "\n";
    return s;
}

void TAntHeader::Print(Option_t*) const
{
    cout << *this;
}

void TAntHeader::Print() const
{
    Print("");
}

