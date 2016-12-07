{
	TH2* data;
	TH2* curr_mc;

	_file0->GetObject("ETheta/sigma", data);
	_file1->GetObject("ETheta/sigma", curr_mc);

	TFile* f = new TFile("step.root","recreate");
	TH2* new_s = ant::TwoPi0_MCSmearing_Tool::CalculateInitialSmearing(data, curr_mc);
	f->Write();

}

