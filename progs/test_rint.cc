
#include "TRint.h"
#include "TCanvas.h"

#include <iostream>

using namespace std;

struct MyCanvas : TCanvas {
    using TCanvas::TCanvas;
};

struct MyWorker  {

    unsigned i = 0;
    TCanvas* c = nullptr;

    void DoWork() {

        if(i<10) {
            i++;
        }
        else if(i==10) {
            c = new TCanvas("bla", "bla");
            c->Draw();
        }
    }

};

struct MyApp : TRint {
    int argc = 0;
    MyWorker worker;

    MyApp(int* argc, char** argv, const MyWorker& worker_) :
        TRint("MyApp", argc, argv, nullptr, 0, false),
        worker(worker_)
    {

    }

    virtual void StartIdleing() override {
        TRint::StartIdleing();
        //cout << "StartIdleing" << endl;
        //worker.DoWork();
    }

};

int main(int argc, char** argv) {


    MyWorker worker;
    MyApp app(&argc, argv, worker);

    app.Run(kTRUE);


}
