#include "TAntHeader.h"

#include <iostream>

#include "TDirectory.h"

#include "base/Logger.h"

using namespace std;
using namespace ant;

TAntHeader::TAntHeader(const std::string& title):
    TNamed("AntHeader", title.c_str()) {
}

bool TAntHeader::IsCompatible(const TAntHeader& other) const {
    return SetupName == other.SetupName
            && GitInfo == other.GitInfo
            && GitInfoDatabase == other.GitInfoDatabase;
}

namespace ant {
ostream& operator<<(ostream& s, const TAntHeader& o) {
    s << "TAntHeader:\n"
      << "  FirstID:         " << o.FirstID << "\n"
      << "  LastID:          " << o.LastID  << "\n"
      << "  Setup:           " << o.SetupName << "\n"
      << "  CmdLine:         " << o.CmdLine << "\n"
      << "  WorkingDir:      " << o.WorkingDir << "\n"
      << "  GitInfo:         " << o.GitInfo << "\n"
      << "  GitInfoDatabase: " << o.GitInfoDatabase << "\n";
    return s;
}
}

void TAntHeader::Print(Option_t*) const
{
    cout << *this;
}

void TAntHeader::Print() const
{
    Print("");
}

void TAntHeader::Browse(TBrowser*)
{
    Print("");
}

Long64_t TAntHeader::Merge(TCollection* li)
{
    if(!li)
        return 0;

    if(FirstID.IsInvalid() || LastID.IsInvalid()) {
        LOG(ERROR) << "Header with invalid FirstID or LastID encountered";
        return 0;
    }

    TIter next(li);

    while(auto h = dynamic_cast<TAntHeader*>(next())) {
        if(!IsCompatible(*h)) {
            LOG(ERROR) << "Skipping incompatible header:\n "
                         << *this << "\n"
                         << *h;
            continue;
        }
        if(h->FirstID.IsInvalid() || h->LastID.IsInvalid()) {
            LOG(ERROR) << "Skipping header with invalid FirstID or LastID:\n"
                       << *h;

            continue;
        }
        FirstID = std::min(FirstID, h->FirstID);
        LastID = std::max(LastID, h->LastID);
    }

    return 0;
}

