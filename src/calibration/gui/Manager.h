#pragma once

#include "AvgBuffer_traits.h"

#include <memory>
#include <list>
#include <vector>
#include <string>

class TH1;
class TFile;
class TQObject;

namespace ant {
namespace calibration {
namespace gui {

class CalCanvasMode;
class ManagerWindowGUI_traits;
class CalibModule_traits;

class Manager {

protected:

    std::unique_ptr<CalibModule_traits>    module;
    std::unique_ptr<AvgBuffer_traits<TH1>> buffer;

    struct input_file_t {

        input_file_t(const std::string& FileName, const interval<TID>& R) :
            filename(FileName),
            range(R)
        {}

        bool operator< (const input_file_t& other) const;

        const std::string filename;
        const ant::interval<ant::TID> range;
    };

    std::list<input_file_t> input_files;

    struct state_t {
        // set in DoInit()
        std::list<input_file_t>::iterator it_file;
        int channel;
        unsigned slice;
        bool oneslice = false;

        bool breakpoint_fit = false;
        bool breakpoint_finish = false;
    };
    state_t state;


    ManagerWindowGUI_traits* window = nullptr;

    void BuildInputFiles(const std::vector<std::string>& filenames);

    void FillBufferFromFiles();

    int nChannels;

    bool confirmed_HeaderMismatch = false;

public:
    std::string SetupName;

    Manager(const std::vector<std::string>& inputfiles,
            std::unique_ptr<AvgBuffer_traits<TH1>> buffer_,
            bool confirmHeaderMismatch=false);

    void SetModule(std::unique_ptr<CalibModule_traits> module_);

    bool DoInit(int gotoSlice);
    void InitGUI(ManagerWindowGUI_traits* window_);

    enum class RunReturn_t {
        Continue, Wait, Exit
    };

    RunReturn_t Run();

    ~Manager();

};

}
}
}
