#pragma once

#include "Manager_traits.h"

#include "base/interval.h"

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
class ManagerWindow;
class AvgBuffer_traits;

class Manager {

protected:

    std::shared_ptr<CalibModule_traits> module;
    std::unique_ptr<AvgBuffer_traits>   buffer;

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
        // set in DoInit()
        std::list<input_file_t>::iterator it_file;
        int channel;
        unsigned slice;
        bool oneslice = false;

        bool breakpoint_fit = false;
        bool breakpoint_finish = false;
    };
    state_t state;


    ManagerWindow* window = nullptr;

    void BuildInputFiles(const std::vector<std::string>& filenames);

    void FillBufferFromFiles();

    int nChannels;

    bool confirmed_HeaderMismatch = false;

public:
    std::string SetupName;

    Manager(const std::vector<std::string>& inputfiles,
            std::unique_ptr<AvgBuffer_traits> buffer_,
            bool confirmHeaderMismatch=false);

    void SetModule(const std::shared_ptr<CalibModule_traits>& module_) {
        module = move(module_);
    }

    bool DoInit(int gotoSlice);
    void InitGUI(ManagerWindow* window_);

    enum class RunReturn_t {
        Continue, Wait, Exit
    };

    RunReturn_t Run();

    ~Manager();

};

}
}
}
