#include "Editor.h"

#include "base/WrapTFile.h"


using namespace std;
using namespace ant;
using namespace ant::calibration;

Editor::Editor(const string& fileName): filename(fileName)
{
    ResetData();
}

void Editor::ResetData()
{
    WrapTFileInput wfi(filename);
    wfi.GetObjectClone("cdata",cdata);
}


void Editor::Save() const
{
    WrapTFileOutput wfo(filename);
    wfo.WriteObject(addressof(cdata),"cdata");
}

void Editor::SaveAs(const string& currentFileName) const
{
    WrapTFileOutput wfo(currentFileName);
    wfo.WriteObject(addressof(cdata),"cdata");
}
