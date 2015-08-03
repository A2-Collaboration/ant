#include "expconfig/ExpConfig.h"
#include "analysis/OutputManager.h"

#include "base/Logger.h"

#include "TRint.h"
#include "TFile.h"
#include "TH1D.h"
#include <iostream>

using namespace std;
using namespace ant;

int main(int argc, char** argv) {
    SetupLogger();
    el::Loggers::setVerboseLevel(9);

    ExpConfig::ManualSetupName = "Setup_2014_EtaPrime";

    auto setup = ExpConfig::Setup::GetLastFound();

    cout << setup->GetName() << endl;

    auto calibs = setup->GetCalibrations();


//    output::OutputManager om;
//    om.SetNewOutput("test.root");

    TRint app("test",&argc,argv,nullptr,0,true);

    WrapTFile file("test.root");

    auto module = calibs.front()->GetPhysicsModule();



//    TH1D* h = new TH1D("h","h",100,0,100);
//    h->Fill(50);



//    h->Draw();

    module->ShowResult();

    app.Run(kTRUE); // really important to return


    cout << "File open=" << file.IsOpen() << endl;

//    file->Write();
//    file->Close();
}
