{
	TH2* data;
	TH2* curr_mc;
	TH2* prev_s;

	_file0->GetObject("ETheta/sigma", data);
	_file1->GetObject("ETheta/sigma", curr_mc);
	_file2->GetObject("energy_smearing", prev_s);

	TFile* f = new TFile("step.root","recreate");
	TH2* new_s = ant::TwoPi0_MCSmearing_Tool::CalculateUpdatedSmearing(data, curr_mc, prev_s);
	f->Write();

}

