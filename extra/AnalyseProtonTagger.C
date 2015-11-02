{
    const char* ofn = Form("%s.root",_file0->GetName());
    cout << ofn << endl;
    
    ProtonTagger->cd();

    TTree* tree = tree;
 
    TFile* of = new TFile(ofn,"RECREATE");

    new TCanvas();
    tree->Draw("tagTime>>h_tagtime(200,-100,100)");

    TCut prompt = "tagTime>-5 && tagTime<5"; // w=10
    TCut random = "(tagTime>-50 && tagTime<-10) || (tagTime>10 && tagTime<50)"; // w = 2 * 40 = 80
    const double ratio = 10.0/80.0;

    //TCut c_ggIM = "ggIM>125 && ggIM<145";
    //TCut c_ggIM = "ggIM>520 && ggIM<570";
    TCut c_ggIM = "ggIM>940 && ggIM<1000";

    TCut cut = c_ggIM;

    new TCanvas();
    tree->Draw("MM>>h_mm_prompt(2000,0,2000)",cut+prompt);

    new TCanvas();
    tree->Draw("MM>>h_mm_random(2000,0,2000)",cut+random);

    new TCanvas();
    TH1F* h_mm = (TH1F*) h_mm_prompt->Clone();
    h_mm->SetName("h_mm");
    h_mm->Add(h_mm_random,-ratio);
    h_mm->Draw();

    TCut c_mm = "MM>800 && MM<1000";

    cut += c_mm;

    new TCanvas();
    tree->Draw("angle>>h_angle_p(200,0,80)",prompt+cut);

    new TCanvas();
    tree->Draw("angle>>h_angle_r(200,0,80)",random+cut);

    new TCanvas();
    TH1F* h_angle = (TH1F*) h_angle_p->Clone();
    h_angle-SetName("h_angle");
    h_angle->Add(h_angle_r,-ratio);
    h_angle->Draw();

    TCut c_angle = "angle<10.0";
    cut += c_angle;

    new TCanvas();
    tree->Draw("E:time>>h_tof_p(100,-20,20,100,0,1000)", prompt+cut,"colz");

    new TCanvas();
    tree->Draw("E:time>>h_tof_r(100,-20,20,100,0,1000)", random+cut,"colz");


    new TCanvas();
    TH2F* h_tof = (TH2F*) h_tof_p->Clone();
    h_tof->SetName("h_tof");
    h_tof->Add(h_tof_r,-ratio);
    h_tof->Draw("colz");

    of->Write();

}

