void plot(const char* pname, const int Cal) {
//	const char* pname = "Proton";
//	const int   Cal   = 2; // 1 = CB, 2 = TAPS, 0 = someone fucked up

	const TCut basecut = Form("mult==1 && rCal == %d", Cal);

	ReconstructCheck->cd();

	new TCanvas();
	tree->Draw("rVeto:rE>>dEE(500,0,500,100,0,10)",basecut,"colz");
	dEE->SetTitle(Form("TAPS %s dEE", pname));
	dEE->GetXaxis()->SetTitle("E [MeV]");
	dEE->GetYaxis()->SetTitle("Veto [MeV]");

	new TCanvas();
	tree->Draw("rSize:rE>>clsz(500,0,500,19,1,20)",basecut,"colz");
	clsz->SetTitle(Form("TAPS %s Cluster Size", pname));
	clsz->GetXaxis()->SetTitle("E [MeV]");
	clsz->GetYaxis()->SetTitle("ClusterSize (Crystals)");

	new TCanvas();
	tree->Draw("rSize:rE>>clsz_he(500,0,500,19,1,20)",basecut + "tE>400","colz");
	clsz_he->SetTitle(Form("TAPS %s Cluster Size (Punching)", pname));
	clsz_he->GetXaxis()->SetTitle("E [MeV]");
	clsz_he->GetYaxis()->SetTitle("ClusterSize (Crystals)");

	new TCanvas();
	tree->Draw("rSize:rE>>clsz_le(500,0,500,19,1,20)", basecut + "tE<400","colz");
	clsz_le->SetTitle(Form("TAPS %s Cluster Size (Non-Punching)", pname));
	clsz_le->GetXaxis()->SetTitle("E [MeV]");
	clsz_le->GetYaxis()->SetTitle("ClusterSize (Crystals)");


	new TCanvas();
	tree->Draw("rE:rTime-5>>tof(400,-20,20,500,0,500)",basecut,"colz");
	tof->SetTitle(Form("TAPS %s ToF", pname));
	tof->GetXaxis()->SetTitle("Time [ns]");
	tof->GetYaxis()->SetTitle("E [MeV]");

}
