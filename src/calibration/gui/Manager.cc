#include "Manager.h"

#include "ManagerWindow.h"

#include "calibration/gui/CalCanvas.h"

#include "tree/TDataRecord.h"
#include "tree/TAntHeader.h"

#include "base/interval.h"
#include "base/std_ext.h"
#include "base/WrapTFile.h"
#include "base/Logger.h"

#include "TH2D.h"

#include <memory>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;


void Manager::BuildInputFiles(const vector<string>& filenames)
{
    if(filenames.empty())
        return;

    for(const auto& filename : filenames) {

        try {

            WrapTFileInput file;
            file.OpenFile(filename);

            TAntHeader* header = nullptr;
            file.GetObject("AntHeader", header);

            if(!header) {
                LOG(WARNING) << "No TAntHeader found in " << filename;
                continue;
            }

            if(SetupName.empty()) {
                SetupName = header->SetupName;
            }
            else if(SetupName != header->SetupName) {
                LOG(WARNING) << "Previously found setup name '" << SetupName
                             << "' does not match '" << header->SetupName << "' of file "
                             << filename;
                continue;
            }

            auto range = interval<TID>(header->FirstID, header->LastID);
            if(!range.IsSane()) {
                LOG(WARNING) << "Range " << range << " not sane in " << filename;
                continue;
            }
            input_files.emplace_back(filename, range);

        } catch (const std::runtime_error& e) {
            LOG(WARNING) << "Can't open " << filename << " " << e.what();
        }
    }

    input_files.sort();
    LOG(INFO) << "Loaded " << input_files.size()
              << " files from " << filenames.size() << " provided filenames ("
              << 100*(double)input_files.size()/filenames.size() << " %)";
}



void Manager::FillBufferFromFiles()
{
    while(buffer.Empty() && state.it_file != input_files.end()) {
        const input_file_t& file_input = *state.it_file;
        try
        {
            WrapTFileInput file;
            file.OpenFile(file_input.filename);

            auto hist = file.GetSharedHist(module->GetHistogramName());

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
        buffer.Finish();
    }
}


Manager::Manager(const std::vector<std::string>& inputfiles, unsigned avglength):
    buffer(avglength),
    state()
{
    BuildInputFiles(inputfiles);
}

void Manager::InitGUI(ManagerWindow* window_) {
    window = window_;
    module->InitCanvases(window);
}

void Manager::ConnectReturnFunc(const char* receiver_class, void* receiver, const char* slot)
{
    signalConnection.receiver_class = receiver_class;
    signalConnection.receiver = receiver;
    signalConnection.slot = slot;
}

bool Manager::input_file_t::operator <(const Manager::input_file_t& o) const {
    return range.Start() < o.range.Start();
}

bool Manager::DoInit()
{
    state.channel = 0;
    state.it_file = input_files.begin();

    maxChannels = module->GetNumberOfChannels();
    if(maxChannels==0) {
        LOG(WARNING) << "Module reports zero channels, nothing to do then.";
        return false;
    }

    FillBufferFromFiles();
    if(buffer.Empty()) {
        LOG(WARNING) << "Did not process anything";
        return false;
    }

    state.is_init = true;
    module->StartRange(buffer.CurrentID());

    return true;
}


bool Manager::Run()
{
    if(!state.is_init) {
       if(!DoInit())
           return false;
    }

    if(!state.breakpoint_finish && state.channel < maxChannels) {
        bool noskip = true;
        if(!state.breakpoint_fit) {
            const string& title = std_ext::formatter() << "Channel=" << state.channel
                                                       << " " << buffer.CurrentID();
            buffer.CurrentSum()->SetTitle(title.c_str());

            const auto ret = module->DoFit(buffer.CurrentSum(), state.channel);
            noskip = ret != Manager_traits::DoFitReturn_t::Skip;

            if(ret == Manager_traits::DoFitReturn_t::Display
               || (!window->Mode.autoContinue && noskip)
               ) {
                VLOG(7) << "Open GUI...";
                module->DisplayFit();
                state.breakpoint_fit = true;
                return false;
            }
        }
        if(noskip)
            module->StoreFit(state.channel);
        state.breakpoint_fit = false;
    }

    if(state.breakpoint_finish
       || (state.channel >= maxChannels && window->Mode.gotoNextSlice))
    {

        if(!state.breakpoint_finish) {
            VLOG(7) << "Finish module first";
            if(module->FinishRange()) {
                VLOG(7) << "GUI Opened (finish range)";
                state.breakpoint_finish = true;
                return false;
            }
        }

        module->StoreFinishRange(buffer.CurrentID());
        state.breakpoint_finish = false;
        state.channel = 0;

        buffer.GotoNextID();

        // try refilling the worklist
        FillBufferFromFiles();
        if(buffer.Empty()) {
            /// \todo give module a chance to do something again here???
            for(CalCanvas* canvas : module->GetCanvases()) {
                canvas->Close();
            }
            LOG(INFO) << "Finished processing whole buffer";
            return false;
        }

        module->StartRange(buffer.CurrentID());

    }
    else
    {
        state.channel += window->Mode.channelStep;
        if(state.channel<0)
            state.channel = 0;
        else if(state.channel>=maxChannels && !window->Mode.gotoNextSlice)
            state.channel = module->GetNumberOfChannels()-1;
    }

    // continue running by default
    return true;
}

Manager::~Manager()
{

}
