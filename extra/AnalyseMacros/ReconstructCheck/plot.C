const char* detFromInt(const int cal) {
    switch (cal) {
        case 1:
            return "CB";
        case 2:
            return "TAPS";
        default:
            return "---";
    }
}
void plot(const char* pname, const int Cal) {
//	const char* pname = "Proton";
//	const int   Cal   = 2; // 1 = CB, 2 = TAPS, 0 = someone fucked up

	const TCut basecut = Form("mult==1 && rCal == %d", Cal);
    const char* det = detFromInt(Cal);

	ReconstructCheck->cd();

	 new TCanvas("cdEE","dEE");
	tree->Draw("rVeto:rE>>dEE(500,0,500,100,0,10)",basecut,"colz");
	dEE->SetTitle(Form("%s %s dEE", det, pname));
	dEE->GetXaxis()->SetTitle("E [MeV]");
	dEE->GetYaxis()->SetTitle("Veto [MeV]");

	new TCanvas("cclsz","Cluster Size");
	tree->Draw("rSize:rE>>clsz(500,0,500,19,1,20)",basecut,"colz");
	clsz->SetTitle(Form("%s %s Cluster Size", det, pname));
	clsz->GetXaxis()->SetTitle("E [MeV]");
	clsz->GetYaxis()->SetTitle("ClusterSize (Crystals)");

	new TCanvas("cclsz_he","Cluster Size (HE)");
	tree->Draw("rSize:rE>>clsz_he(500,0,500,19,1,20)",basecut + "tE>400","colz");
	clsz_he->SetTitle(Form("%s %s Cluster Size (Punching)", det, pname));
	clsz_he->GetXaxis()->SetTitle("E [MeV]");
	clsz_he->GetYaxis()->SetTitle("ClusterSize (Crystals)");

	new TCanvas("cclsz_le","Cluster Size (LE)");
	tree->Draw("rSize:rE>>clsz_le(500,0,500,19,1,20)", basecut + "tE<400","colz");
	clsz_le->SetTitle(Form("%s %s Cluster Size (Non-Punching)", det, pname));
	clsz_le->GetXaxis()->SetTitle("E [MeV]");
	clsz_le->GetYaxis()->SetTitle("ClusterSize (Crystals)");

	new TCanvas("ctof","ToF");
	tree->Draw("rE:rTime-5>>tof(400,-20,20,500,0,500)",basecut,"colz");
	tof->SetTitle(Form("%s %s ToF", det, pname));
	tof->GetXaxis()->SetTitle("Time [ns]");
	tof->GetYaxis()->SetTitle("E [MeV]");

	new TCanvas("ceRec","Energy Reconstruction");
	tree->Draw("rE:tE>>eRec(1600,0,1600,1600,0,1600)",basecut,"colz");
	eRec->SetTitle(Form("%s %s Energy Reconstruction", det, pname));
	eRec->GetXaxis()->SetTitle("E_{true} [MeV]");
	eRec->GetYaxis()->SetTitle("E_{rec} [MeV]");
}
