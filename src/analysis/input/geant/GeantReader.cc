#include "GeantReader.h"


#include "base/WrapTFile.h"
#include "base/Logger.h"

#include "TTree.h"


using namespace ant;
using namespace ant::analysis;
using namespace ant::analysis::input;
using namespace ant::analysis::data;
using namespace std;

GeantReader::GeantReader(const std::shared_ptr<WrapTFileInput>& rootfiles) :
    files(rootfiles)
{

    if(!files->GetObject("h12", tree))
        return;

    VLOG(5) << "Found Geant tree h12";

    tree->SetBranchAddress("vertex",fvertex);

    LOG(INFO) << "MCTrue from geant input active";
}

GeantReader::~GeantReader() {}

bool GeantReader::ReadNextEvent(Event& event)
{
    if(!tree)
        return false;

    if(current_entry >= tree->GetEntries())
        return false;

    tree->GetEntry(current_entry);

    event.MCTrue.Target.Vertex = TVector3(fvertex[0], fvertex[1], fvertex[2]);

    ++current_entry;

    return true;
}

double GeantReader::PercentDone() const
{
    return double(current_entry) / double(tree->GetEntries());
}
