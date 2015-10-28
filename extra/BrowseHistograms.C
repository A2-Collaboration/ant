// little histogram browser,
// useful for searching jumps



const char fFileDir[] = "analysis";
//const char fHistName[] = "CaLib_CB_Time_Ind";
const char fHistName[] = "CaLib_TAPS_Time_Ind";
//const char fHistName[] = "CaLib_Veto_Time_Ind";
//const char fHistName[] = "CaLib_PID_Time_Ind";
//const char fHistName[] = "CaLib_Tagger_Time_Ind";

TList* histograms = 0;
TList* saved = 0;
Int_t currPos = 0;
TCanvas* c = 0;
TFile* currFile = 0;
TH1* currHist = 0;

void NextAndForget() {
  Int_t i_found = -1;
  for(Int_t i=0;i<saved->GetEntries(); i++) {
    if(currHist->GetTitle() == ((TObjString*)saved->At(i))->GetString()) {
      i_found = i;
      break;
    }
  }
  if(i_found<0) {
    cout << "Not found in saved list" << endl;
  }
  else {
    saved->RemoveAt(i_found);
    cout << "Removed from list" << endl;
  }
  ProcessClick(1);
}

void NextAndSave() {
  saved->Add(new TObjString(currHist->GetTitle()));
  cout << "Saved to list" << endl;
  ProcessClick(1);
}

void PrintSaved() {
  TIter next(saved);
  TObjString* title;
  cout << endl << "Saved histograms: " << saved->GetEntries() << endl << endl;
  while ((title=(TObjString*)next())) {
    cout << title->GetString() << endl;
  }
}

void ProcessClick(Int_t direction)
{
  if(currPos==0 && direction<0)
    return;
  if(currPos==histograms->GetEntries()-1 && direction>0)
    return;
  currPos += direction;
  DrawCurrHistogram();
}

void PreserveAxis(TAxis* axis1, TAxis* axis2) {
  Int_t binmin = axis1->GetFirst();
  Int_t binmax = axis1->GetLast();
  Float_t xmin = axis1->GetBinLowEdge(binmin);
  Float_t xmax = axis1->GetBinLowEdge(binmax);
  Int_t newmin = axis2->FindBin(xmin);
  Int_t newmax = axis2->FindBin(xmax);
  axis2->SetRange(newmin,newmax);
}

void DrawCurrHistogram() {


  TObjString* filename = histograms->At(currPos);
  TString f(fFileDir);
  f.Append("/");
  f.Append(filename->GetString());

  TFile* file = new TFile(f, "read");
  // check file
  if (!file) return;
  if (file->IsZombie()) return;

  TH2D* hist = (TH2D*) gFile->Get(fHistName);

  if (!hist) continue;
  if (!hist->GetEntries()) continue;

  // copy the histogram to the root memory
  //gROOT->cd();
  //TH2D* hist_copy = new TH2D(*hist);
  TString t(filename->GetString());
  t.Append(": ");
  t.Append(hist->GetTitle());
  hist->SetTitle(t);
  //histograms->Add(hist_copy);
  //TH1* h = histograms->At(currPos);
  c->cd(1);
  //gPad->cd();

  hist->Draw("colz");
  // file can be closed now...
  //file->Close();

  if(currHist != 0) {
    PreserveAxis(currHist->GetXaxis(), hist->GetXaxis());
    PreserveAxis(currHist->GetYaxis(), hist->GetYaxis());
  }

  if(currFile != 0)
    currFile->Close();
  currFile = file;
  currHist = hist;
  //gPad->Modified();
  //gPad->Update();

  //cout << "Added histogram from " << filename << endl;

}

void BrowseHistograms() {

  TSystemDirectory dir(fFileDir,fFileDir);
  TList *files = dir.GetListOfFiles();
  if (!files) {
    cerr << "Error: No files found in " << fFileDir << endl;
    return;
  }
  files->Sort();


  histograms = new TList();
  saved = new TList();

  c = new TCanvas("canvas","Browse Histograms");

  TSystemFile *file;
  TIter next(files);
  Int_t n=0;
  TString fname;
  while ((file=(TSystemFile*)next())) {
    fname = file->GetName();
    if (!file->IsDirectory() && fname.EndsWith(".root")) {

      histograms->Add(new TObjString(fname));
      n++;
      //if(n>10)
      //	break;
    }
  }

  DrawCurrHistogram();

  TControlBar *bar = new TControlBar("vertical", "Control", 10, 10);
  bar->AddButton("Next","ProcessClick(1);", "Show next run");
  bar->AddButton("Next And Save","NextAndSave();", "Show previous run");
  bar->AddButton("Next And Forget","NextAndForget();", "Show previous run");
  bar->AddButton("Prev","ProcessClick(-1);", "Show previous run");
  bar->AddButton("Print saved","PrintSaved();", "Show previous run");
  bar->SetButtonWidth(90);
  bar->Show();
}
