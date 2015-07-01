#ifndef GOATREADER_H
#define GOATREADER_H

#include "DataReader.h"
#include "Event.h"

#include <memory>
#include <string>
#include <list>

#include "FileManager.h"
#include "TreeManager.h"
#include "InputModule.h"

#include "TriggerInput.h"
#include "TaggerInput.h"
#include "DetectorHitInput.h"
#include "TrackInput.h"

namespace ant {
namespace input {




class GoatReader: public DataReader {
protected:

    class ModuleManager: public std::list<BaseInputModule*> {
    public:
        ModuleManager() = default;
        ModuleManager(const std::initializer_list<BaseInputModule*>& initlist):
            std::list<BaseInputModule*>(initlist) {}

        void GetEntry() {
            for(auto& module : *this) {
                module->GetEntry();
            }
        }

    };

    FileManager   files;
    TreeManager   trees;


    TriggerInput trigger;
    TaggerInput  tagger;

    ModuleManager active_modules = {&trigger, &tagger};

    Long64_t    current_entry = -1;

    void AddInputModule(BaseInputModule& module);

public:
    GoatReader() = default;
    virtual ~GoatReader() = default;

    void AddInputFile(const std::string& filename);
    void Initialize();

    Long64_t  GetNEvents() const;

    std::shared_ptr<Event> ReadNextEvent();
    bool hasData() const;


};

}
}

#endif
