#include "TCanvas.h"
#include "TRint.h"

#include "iostream"
//using namespace ant;
using namespace std;


#include "TGClient.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TRandom.h"
#include "TGButton.h"
#include "TGFrame.h"
#include "TRootEmbeddedCanvas.h"
#include "RQ_OBJECT.h"
#include "KeySymbols.h"

class ExtMainFrame: public TGMainFrame
{
public:
   ExtMainFrame(const TGWindow *p,UInt_t w,UInt_t h):
       TGMainFrame(p,w,h)
   {
       BindKey(p,gVirtualX->KeysymToKeycode(kKey_0),0);
   }
   virtual Bool_t HandleKey(Event_t *event) override
   {
       cout << "pressed  key" << endl;

       char tmp[2];
       UInt_t keysym;

       gVirtualX->LookupString(event, tmp, sizeof(tmp), keysym);



       if (event->fType == kGKeyPress) {
           switch ((EKeySym)keysym) {
              case kKey_Left:
                 cout << "pressed Up key" << endl;
                 break;
              default:
                 break;
           }
        }
       return TGMainFrame::HandleKey(event);
   }








};

class MyMainFrame 
{
   RQ_OBJECT("MyMainFrame")
private:
   ExtMainFrame         *fMain;
//   TRootEmbeddedCanvas *fEcanvas;
public:
   MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
   virtual ~MyMainFrame();
   void DoDraw();
};

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h) 
{
   // Create a main frame
   fMain = new ExtMainFrame(p,w,h);

   // Create canvas widget
   //fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,200,200);
   //fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX |
                   //kLHintsExpandY, 10,10,10,1));
   // Create a horizontal frame widget with buttons
   TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);
//   TGTextButton *draw = new TGTextButton(hframe,"&Draw");
//   draw->Connect("Clicked()","MyMainFrame",this,"DoDraw()");
//   hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX,
//                                            5,5,3,4));
//   TGTextButton *exit = new TGTextButton(hframe,"&Exit",
//                                "gApplication->Terminate(0)");
//   hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,
//                                            5,5,3,4));
//   fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,
//                                             2,2,2,2));

   // Set a name to the main frame
   fMain->SetWindowName("Simple Example");

   // Map all subwindows of main frame
   fMain->MapSubwindows();

   // Initialize the layout algorithm
   fMain->Resize(fMain->GetDefaultSize());

   // Map main frame
   fMain->MapWindow();
//   DoDraw();
}

void MyMainFrame::DoDraw() 
{
   // Draws function graphics in randomly choosen interval
//   TF1 *f1 = new TF1("f1","sin(x)/x",0,gRandom->Rndm()*10);
//   f1->SetLineWidth(3);
//   f1->Draw();
//   TCanvas *fCanvas = fEcanvas->GetCanvas();
//   fCanvas->cd();
//   fCanvas->Update();
}

MyMainFrame::~MyMainFrame() 
{
   // Clean up used widgets: frames, buttons, layout hints
   fMain->Cleanup();
   delete fMain;
}

int main(int argc, char** argv) 
{
   int i = 0;
   TRint app("guitest",&i,nullptr);
   new MyMainFrame(gClient->GetRoot(),200,200);
   app.Run(kFALSE);
}
