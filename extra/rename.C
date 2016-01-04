{
    _file0->ReOpen("UPDATE");
    
    TDirectory* d = NULL;

    gDirectory->GetObject("TAPS_Energy", d);

    if(!d) {
        cout << "TAPS_Energy not found";
        exit(2);
    }

    d->cd();

    TH2* h_mm = NULL;
    h_mm = d->GetObject("ggIM",h_mm);

    if(h_mm) {
        h_mm->SetName("RelativeGains");
    }
    else
        cout << "ggIM not found" << endl;

    _file0->Write();

   exit(0);

}
