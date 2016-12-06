//
// Fit pi0 peak in each bin
// root <root-file> Analyse.C
//  reates root-file.a.root with histograms
//
{
	TH3* h_ETheta;
	_file0->GetObject("TwoPi0_MCSmearing/m2Pi0/cb_pi0_ETheta", h_ETheta);
	if(!h_ETheta) {
		cerr << "Hist not found!" << endl;
	}

//	TH3* h_EElement;
//	_file0->GetObject("TwoPi0_MCSmearing/m2Pi0/cb_pi0_E_Element", h_EElement);
//	if(!h_EElement) {
//		cerr << "Hist not found!" << endl;
//	}

	
//	TH2* h_Element;
//	_file0->GetObject("TwoPi0_MCSmearing/m2Pi0/cb_pi0", h_Element);
//	if(!h_Element) {
//		cerr << "Hist not found!" << endl;
//	}

	TFile* f = new TFile(Form("%s.a.root", _file0->GetName()), "recreate");

	TDirectory* base = gDirectory;
	
//	TDirectory* d = base->mkdir("Element");
//	d->cd();
//	ant::TwoPi0_MCSmearing_Tool::FitAllChannels(h_Element);

	d = base->mkdir("ETheta");
	d->cd();
	ant::TwoPi0_MCSmearing_Tool::AnalyseChannelE(h_ETheta);

//	d = base->mkdir("EElement");
//	d->cd();
//	ant::TwoPi0_MCSmearing_Tool::AnalyseChannelE(h_EElement);
	base->cd();

	f->Write();
}
