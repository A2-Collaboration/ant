#include "CalibrationGUI.h"

#include "AvgBuffer.h"
#include "GUIInterface.h"

#include "calibration/gui/CalCanvas.h"

#include "tree/TDataRecord.h"
#include "tree/TAntHeader.h"

#include "base/interval.h"
#include "base/std_ext.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"

#include "TH1D.h"
#include "TH2D.h"

#include <memory>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;


std::list<CalibrationGUI::input_file_t> CalibrationGUI::ScanFiles(const std::vector<string> filenames)
{
    std::list<CalibrationGUI::input_file_t> inputs;

    for(auto& filename : filenames) {

        try {

            WrapTFile file(filename, WrapTFile::mode_t::read, false);

            TAntHeader* header = nullptr;
            file.GetObject("AntHeader", header);


            if(header) {
                auto i = ant::interval<TID>(header->FirstID, header->LastID);
                LOG(WARNING) << i;
                inputs.emplace_back(filename, i);
            } else {
                LOG(WARNING) << "No TAntHeader found in " << filename;
            }
        } catch (const std::runtime_error& e) {
            LOG(WARNING) << "Can't open " << filename << " " << e.what();

        }
    }

    return inputs;
}



void CalibrationGUI::FillWorklistFromFiles()
{
    while(buffer.Worklist().empty() && state.it_file != input_files.end()) {
        const input_file_t& file_input = *state.it_file;
        try
        {
            WrapTFile file(file_input.filename, WrapTFile::mode_t::read, false);

            auto hist = file.GetSharedTH2(module->GetHistogramName());

            if(!hist) {
                LOG(WARNING) << "Histogram " << module->GetHistogramName() << " not found in " << file_input.filename;
            } else {
                buffer.Push(hist, file_input.range);
            }

        }
        catch (const std::runtime_error& e) {
            LOG(WARNING) << "Can't open " << file_input.filename << ": " << e.what();
        }
        state.it_file++;
    }

    if(state.it_file == input_files.end()) {
        VLOG(7) << "Reached end of files, processing remaining buffer";
        buffer.PushRestToWorklist();
    }
}

CalibrationGUI::CalibrationGUI(std::unique_ptr<GUIClientInterface> module_, unsigned length):
    module(move(module_)),
    buffer(length),
    state(),
    mode(std_ext::make_unique<CalCanvasMode>())
{
    module->InitGUI();
    for(CalCanvas* canvas : module->GetCanvases()) {
        canvas->LinkGUIMode(mode.get());
    }
}

void CalibrationGUI::ConnectReturnFunc(const char* receiver_class, void* receiver, const char* slot)
{
    for(CalCanvas* canvas : module->GetCanvases()) {
        canvas->ConnectReturnFunc(receiver_class, receiver, slot);
    }
}

void CalibrationGUI::SetFileList(const std::vector<string>& filelist)
{
    VLOG(7) << "Scanning input files...";
    input_files = ScanFiles(filelist);
    VLOG(7) << "Sorting input files by TID range...";
    input_files.sort();
    VLOG(7) << "Input files scanned";
}

bool CalibrationGUI::input_file_t::operator <(const CalibrationGUI::input_file_t& o) const {
    return range.Start() < o.range.Start();
}

bool CalibrationGUI::Run()
{
    if(!state.is_init) {
        state.channel = 0;
        state.it_buffer = buffer.begin();
        state.it_file = input_files.begin();
        FillWorklistFromFiles();
        if(buffer.Worklist().empty()) {
            LOG(WARNING) << "Did not process anything";
            return false;
        }
        state.is_init = true;
    }

    VLOG(8) << "Worklist size " << buffer.Worklist().size();

    if(!state.stop_finish) {
        if(!state.stop_fit) {
            const string& title = std_ext::formatter() << "Channel=" << state.channel
                                                       << " " << buffer.Worklist().front();
            buffer.Average()->SetTitle(title.c_str());

            const bool stop = module->Fit(buffer.Average(), state.channel);

            if(stop || mode->stopAlways) {
                VLOG(7) << "Open GUI...";
                module->DisplayFit();
                state.stop_fit = true;
                return false;
            }
        }
        module->StoreResult(state.channel);
        state.stop_fit = false;
    }

    if(
       (state.channel >= (int)module->GetNumberOfChannels() && mode->gotoNextBuffer)
       || state.stop_finish
       )
    {

        if(!state.stop_finish) {
            VLOG(7) << "Finish module first";
            if(module->Finish()) {
                VLOG(7) << "GUI Opened (finish)";
                state.stop_finish = true;
                return false;
            }
        }

        module->StoreFinish();
        state.stop_finish = false;
        state.channel = 0;

        buffer.Worklist().pop();

        // try refilling the worklist
        FillWorklistFromFiles();
        if(buffer.Worklist().empty()) {
            /// \todo give module a chance to do something again here???
            for(CalCanvas* canvas : module->GetCanvases()) {
                canvas->Close();
            }
            LOG(INFO) << "Finished processing whole buffer";
            return false;
        }

    }
    else
    {
        state.channel += mode->channelStep;
        if(state.channel<0)
            state.channel = 0;
        else if(state.channel>(int)module->GetNumberOfChannels())
            state.channel = module->GetNumberOfChannels();
    }

    // continue running by default
    return true;
}

CalibrationGUI::~CalibrationGUI()
{

}
