
#include "base/std_ext.h"

#include "TRint.h"
#include "TCanvas.h"
#include "TThread.h"
#include "RQ_OBJECT.h"
#include "TExec.h"

#include <iostream>

using namespace std;
using namespace ant;


// should be replaced by real ant gui worker class
struct MyWorker {
    unsigned i = 0;
    TCanvas* DoWork() {
        for(;i<100;i++) {
            cout << "Work i=" << i << endl;

            if(i==10 || i==20) {
                cout << "Canvas needed i=" << i << endl;
                i++;
                return new TCanvas("bla", "bla");
            }
        }
        cout << "DONE" << endl;
        return nullptr;
    }
};


struct MyExec : TExec {

    //TCanvas* c = nullptr;
    std::unique_ptr<MyWorker> worker;

    MyExec(std::unique_ptr<MyWorker> worker_) :
        worker(move(worker_))
    {}

    virtual void Exec(const char*) override {
        TCanvas* c = worker->DoWork();
        if(c == nullptr)
            return;
        c->Connect("Destroyed()", "TExec", this, "Exec(=\"\")");
        c->Draw();
    }
};

int main(int argc, char** argv) {

    auto app = std_ext::make_unique<TRint>("app",&argc,argv);
    auto exec = std_ext::make_unique<MyExec>(std_ext::make_unique<MyWorker>());
    exec->Exec("");
    app->Run(kTRUE);
    app = nullptr;
    exec = nullptr;
}
