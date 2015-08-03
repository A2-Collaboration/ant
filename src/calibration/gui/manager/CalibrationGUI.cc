#include "CalibrationGUI.h"
#include "AvgBuffer.h"
#include "base/interval.h"
#include "tree/TDataRecord.h"
#include "TH1.h"
#include "TH2.h"
#include "base/std_ext.h"
#include "TFile.h"
#include "base/Logger.h"

#include <memory>

using namespace std;
using namespace ant;
using namespace ant::calibration;
using namespace ant::calibration::gui;

class CalibrationGUI::myBuffer_t: public AvgBuffer<TH1,ant::interval<ant::TID>> {
    using AvgBuffer<TH1,ant::interval<ant::TID>>::AvgBuffer;
};

void CalibrationGUI::ReadFile(const std::string& filename)
{
    auto file = std_ext::make_unique<TFile>(filename.c_str(),"READ");

    //@todo: ant histgoram file reader

    if(!file->IsZombie()) {

        TID* first_it = nullptr;
        file->GetObject("", first_it);

        TID* last_it = nullptr;
        file->GetObject("", last_it);

        if(last_it && first_it) {

            for(auto& mod : mod_buffers) {

                TH1* hist = nullptr;
                file->GetObject(mod.module->GetHistogramName().c_str(), hist);

                if(!hist) {
                    LOG(WARNING) << "Histogram " << mod.module->GetHistogramName() << " not found in " << filename;
                } else {
                   auto shist = shared_ptr<TH1>(static_cast<TH1*>(hist->Clone()));
                   mod.buffer->Push(shist, interval<TID>(*first_it, *last_it));
                }
            }
        } else {
            LOG(WARNING) << "TID interval not found/complete in " << filename;
        }
    }

    file->Close();
}

void CalibrationGUI::ProcessModules()
{
    for(auto& mod : mod_buffers) {
        if(mod.buffer->isFull()) {
            // process histograms
        }
    }
}

void CalibrationGUI::Run(const std::vector<std::string>& filelist)
{
    for(auto& filename : filelist) {
        ReadFile(filename);
        ProcessModules();
    }
}

void CalibrationGUI::AddModule(GUIClientInrerface* module, unsigned avg_length)
{
    mod_buffers.emplace_back( module, move(std_ext::make_unique<myBuffer_t>(avg_length)));
}

CalibrationGUI::~CalibrationGUI()
{

}
