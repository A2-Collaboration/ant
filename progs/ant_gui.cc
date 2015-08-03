
#include "base/std_ext.h"

#include "TRint.h"
#include "TCanvas.h"
#include "TThread.h"
#include "RQ_OBJECT.h"
#include "TExec.h"
#include "TObject.h"
#include "base/Logger.h"
#include "calibration/gui/debug_GUIModule.h"
#include "calibration/gui/manager/CalibrationGUI.h"
#include "base/CmdLine.h"

#include <iostream>

using namespace std;
using namespace ant;


//// should be replaced by real ant gui worker class
//struct MyWorker {
//    struct state_t {
//        int i;
//        bool i_init;
//        int j;
//        bool break_occured;
//    };

//    state_t t = {0,false,0,false};
//    MyWorker() {
//        t.i = 0;
//        t.i_init = false;
//        t.j = 0;
//        t.break_occured = false;
//    }

//    TCanvas* DoWork() {

//        while(t.i <100) {

//            if(!t.i_init) {
//                t.i_init = true;
//                cout << "INIT i" << t.i << endl;
//            }

//            while(t.j<100) {

//                if(t.break_occured) {
//                    cout << "returning from break. j="<<t.j << endl;
//                    t.break_occured = false;
//                    t.j++;
//                }

//                cout << "i=" << t.i << " " << t.j << endl;

//                if((t.j%10)==0) {
//                    cout << "Need break" <<endl;
//                    t.break_occured = true;
//                    return new TCanvas("bla", "bla");
//                }
//                ++t.j;
//            }
//            t.j=0;
//            t.i_init=false;
//            ++t.i;
//        }

//        cout << "DONE" << endl;
//        return nullptr;
//    }
//};

using namespace  ant::calibration::gui;

struct MyExec : TExec {

    const TRint* rint = nullptr;
    std::unique_ptr<ant::calibration::gui::CalibrationGUI> gui;

    MyExec(std::unique_ptr<ant::calibration::gui::CalibrationGUI> gui_, const TRint* Rint) :
        rint(Rint),
        gui(move(gui_))
    {}

    virtual void Exec(const char* arg) override {

        const string Arg(arg);

        if(Arg != "firstcall" && !rint->IsRunning()) {
            return;
        }

        CalibrationGUI::RunReturn_t c;

        do {

            c = gui->Run();
        } while (c.status == CalibrationGUI::RunReturnStatus_t::Next);

        if( c.status == CalibrationGUI::RunReturnStatus_t::OpenGUI ) {
            c.gui->Connect("Destroyed()", "TExec", this, "Exec(=\"\")");
        }
    }
};

int main(int argc, char** argv) {
    SetupLogger();


    TCLAP::CmdLine cmd("ant", ' ', "0.1");
    auto cmd_verbose = cmd.add<TCLAP::ValueArg<int>>("v","verbose","Verbosity level (0..9)", false, 0,"level");
    auto cmd_input  = cmd.add<TCLAP::MultiArg<string>>("i","input","Input files",true,"filename");
    cmd.parse(argc, argv);

    if(cmd_verbose->isSet())
    el::Loggers::setVerboseLevel(cmd_verbose->getValue());

    ant::calibration::gui::DebugModule d;

    int i=0;
    auto app = std_ext::make_unique<TRint>("app",&i,nullptr);
    std::unique_ptr<ant::calibration::gui::CalibrationGUI> gui = std_ext::make_unique<ant::calibration::gui::CalibrationGUI>(&d,5);
    gui->SetFileList(cmd_input->getValue());
    gui->Prepare();

    auto exec = std_ext::make_unique<MyExec>(std::move(gui), app.get());
    exec->Exec("firstcall");
    app->Run(kTRUE);
    app = nullptr;
    exec = nullptr;
}
