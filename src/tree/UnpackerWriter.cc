#include "UnpackerWriter.h"

#include "TEvent.h"
#include "TDetectorRead.h"
#include "THeaderInfo.h"
#include "TUnpackerMessage.h"
#include "TSlowControl.h"


#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "base/std_ext.h"

#include <algorithm>

using namespace std;
using namespace ant;
using namespace ant::tree;


UnpackerWriter::UnpackerWriter(const string& outputfile)
{
    file = std_ext::make_unique<WrapTFile>(outputfile);
    Event.Init();
    DetectorRead.Init();
    HeaderInfo.Init();
    UnpackerMessage.Init();
    SlowControl.Init();
}

UnpackerWriter::~UnpackerWriter() {}

void UnpackerWriter::Fill(TDataRecord* record) noexcept
{
   if(Event.TryFill(record))
       return;
   if(DetectorRead.TryFill(record))
       return;
   if(HeaderInfo.TryFill(record))
       return;
   if(UnpackerMessage.TryFill(record))
       return;
   if(SlowControl.TryFill(record))
       return;
}
