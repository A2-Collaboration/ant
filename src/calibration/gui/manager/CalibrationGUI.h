#pragma once

#include "GUIInterface.h"
#include <memory>
#include <list>
#include <vector>

class TH1;

namespace ant {
namespace calibration {
namespace gui {

class CalibrationGUI {
protected:

    class myBuffer_t;

    struct module_buffer {
        module_buffer(GUIClientInrerface* gi, std::unique_ptr<myBuffer_t> bptr): module(gi), buffer(move(bptr)) {}
        GUIClientInrerface* module;
        std::unique_ptr<myBuffer_t> buffer;
    };

    std::list<module_buffer> mod_buffers;

    void ReadFile(const std::string& filename);
    void ProcessModules();

public:

    virtual void Run(const std::vector<std::string>& filelist);

    virtual void AddModule(GUIClientInrerface* module, unsigned avg_length);

    virtual ~CalibrationGUI();

};

}
}
}
