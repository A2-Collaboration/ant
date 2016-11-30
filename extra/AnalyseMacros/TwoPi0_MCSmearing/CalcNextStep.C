{
	TH2* data;
	TH2* curr_mc;
	TH2* prev_diff;
	TH2* diff_sum;

	_file0->GetObject("EElement/sigma", data);
	_file1->GetObject("EElement/sigma", curr_mc);
	_file2->GetObject("im_d", prev_diff);
	_file2->GetObject("im_d_sum", diff_sum);

	TFile* f = new TFile("step.root","recreate");
	TH2* new_s = ant::TwoPi0_MCSmearing_Tool::CalculateUpdatedSmearing(data, curr_mc, prev_diff, diff_sum);
	f->Write();

}

