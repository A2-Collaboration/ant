#pragma once

#include "GUIInterface.h"
#include "AvgBuffer.h"

#include "base/interval.h"
#include "tree/TDataRecord.h"

#include <memory>
#include <list>
#include <vector>
#include <string>

class TH1D;
class TH2D;
class TFile;
class TQObject;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvasMode;

class CalibrationGUI {
protected:

    using myBuffer_t = AvgBuffer<TH2D, interval<TID>>;

    std::unique_ptr<GUIClientInterface> module;
    myBuffer_t buffer;

    struct input_file_t {

        input_file_t(const std::string& FileName, const interval<TID>& R):
            filename(FileName),
            range(R) {}

        bool operator< (const input_file_t& other) const;

        const std::string filename;
        const ant::interval<ant::TID> range;
    };

    std::list<input_file_t> input_files;

    struct state_t {
        state_t() :
              is_init(false),
              stop_fit(false),
              stop_finish(false)
        {}

        std::list<input_file_t>::iterator it_file;
        myBuffer_t::const_iterator it_buffer;
        int channel;

        bool is_init;
        bool stop_fit;
        bool stop_finish;
    };
    state_t state;


    std::unique_ptr<CalCanvasMode> mode;

    std::list<input_file_t> ScanFiles(const std::vector<std::string> filenames);

    void FillWorklistFromFiles();

public:
    enum class RunReturnStatus_t {
        Continue,
        Stop
    };

    CalibrationGUI(std::unique_ptr<GUIClientInterface> module_,
                   unsigned length);

    virtual void ConnectReturnFunc(const char* receiver_class, void* receiver, const char* slot);

    virtual void SetFileList(const std::vector<std::string>& filelist);

    virtual bool Run();

    virtual ~CalibrationGUI();

};

}
}
}
