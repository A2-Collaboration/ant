#pragma once

#include <string>
#include <list>

class TH1;
class TQObject;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvas;

class GUIClientInterface {
public:

    virtual std::string GetHistogramName() const =0;
    virtual unsigned GetNumberOfChannels() const =0;
    virtual std::list<CalCanvas*> GetCanvases() const =0;
    virtual void InitGUI() =0;

    virtual bool Fit(TH1* hist, unsigned channel) =0;
    virtual void DisplayFit() =0;
    virtual void StoreResult(unsigned channel) =0;

    virtual bool Finish() =0;
    virtual void StoreFinish() =0;

};

}
}
}
