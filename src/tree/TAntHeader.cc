#include "TAntHeader.h"

#include "TDirectory.h"

ant::TAntHeader::TAntHeader(const std::string& title):
    TNamed("AntHeader", title.c_str()) {
    gDirectory->Add(this);
}

ant::TAntHeader::~TAntHeader()
{
}
