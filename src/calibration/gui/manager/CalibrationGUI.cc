#include "CalibrationGUI.h"
#include "AvgBuffer.h"
#include "base/interval.h"
#include "tree/TDataRecord.h"
#include "TH1D.h"
#include "TH2D.h"
#include "base/std_ext.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"
#include "tree/TAntHeader.h"
#include "gui/FitCanvas.h"
#include "GUIInterface.h"
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
                LOG(WARNING) << "no TAntHeader in " << filename;
            }
        } catch (const std::runtime_error& e) {
            LOG(WARNING) << "Can't open " << filename << " " << e.what();

        }
    }

    return inputs;
}



void CalibrationGUI::ProcessFile(input_file_t& file_input)
{

    try {
        WrapTFile file(file_input.filename, WrapTFile::mode_t::read, false);


        auto hist = file.GetSharedTH2(module->GetHistogramName());

        if(!hist) {
            LOG(WARNING) << "Histogram " << module->GetHistogramName() << " not found in " << file_input.filename;
        } else {
            buffer.Push(hist, file_input.range);
        }

    } catch (const std::runtime_error& e) {
        LOG(WARNING) << "Can't open " << file_input.filename << " " << e.what();

    }
}

CalibrationGUI::CalibrationGUI(std::unique_ptr<GUIClientInterface> module_, unsigned length):
    module(move(module_)),
    canvas(std_ext::make_unique<CalCanvas>("ModuleName")), /// \todo obtain name from module
    buffer(length)
{
    state.is_init = false;
    state.finish_mode = false;
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

CalibrationGUI::RunReturn_t CalibrationGUI::Run()
{
    if(!state.is_init) {
        state.break_occured=false;
        state.channel=0;
        state.buffpos=buffer.begin();
        state.file=input_files.begin();
        ProcessFile(*state.file);

        state.is_init = true;
    }

    if(state.break_occured == true) {
        VLOG(7) << "Returning from GUI";
        module->StoreResult(state.channel);
        state.break_occured = false;
    } else {

        if(!buffer.Worklist().empty()) {
            const string& title = std_ext::formatter() << "Channel=" << state.channel << " " << buffer.Worklist().top();
            buffer.Average()->SetTitle(title.c_str());
            canvas->Clear();
            GUIClientInterface::FitStatus r = module->Fit(canvas.get(), buffer.Average(), state.channel);

            if(r == GUIClientInterface::FitStatus::GUIWait) {
                VLOG(7) << "GUI Opened";
                state.break_occured = true;
                return RunReturn_t(RunReturnStatus_t::OpenGUI, canvas.get());
            }

            module->StoreResult(state.channel);
        }
    }

    VLOG(8) << "Buffer length " << buffer.Worklist().size();

    if(state.channel >= module->GetNumberOfChannels() || buffer.Worklist().empty()) {
        state.channel = 0;

        if(!buffer.Worklist().empty())
            buffer.Worklist().pop();

        if(buffer.Worklist().empty()) {
            if(state.finish_mode) {
                return RunReturn_t(RunReturnStatus_t::Done);
            }

            ++state.file;

            if(state.file == input_files.end()) {

                if(!state.finish_mode) {
                    buffer.PushRestToWorklist();
                    state.finish_mode = true;
                }

            } else {
                ProcessFile(*state.file);
            }
        }
    } else {
        ++state.channel;
    }

    return RunReturn_t(RunReturnStatus_t::Next);
}

CalibrationGUI::~CalibrationGUI()
{

}
