void CopyBack(const TF1* sum, TF1* sig, TF1* bg) {
	bg->SetParameters(sum->GetParameters());
	bg->SetParErrors(sum->GetParErrors());
	sig->SetParameters(&(sum->GetParameters()[bg->GetNpar()]));
	sig->SetParErrors(&(sum->GetParErrors()[bg->GetNpar()]));
}

void SumSetNames(TF1* sum, const TF1* sig, const TF1* bg) {
	int n=0;
	for(int i=0; i<bg->GetNpar(); ++i) {
		sum->SetParName(n++, bg->GetParName(i));
	}
	
	for(int i=0; i<sig->GetNpar(); ++i) {
		sum->SetParName(n++, sig->GetParName(i));
	}

}

void CopyParLimits(TF1* sum, const TF1* sig, const TF1* bg) {
	
	int n=0;
	double a,b;

	for(int i=0; i<bg->GetNpar(); ++i) {
		bg->GetParLimits(i, a, b);
		const double v = bg->GetParameter(i);

			if(a==0.0 && b==0.0) {
				sum->ReleaseParameter(n);
			} else {
				sum->SetParLimits(n, a, b);
			}
			++n;
	}
	
	for(int i=0; i<sig->GetNpar(); ++i) {
		sig->GetParLimits(i, a, b);
		const double v = sig->GetParameter(i);

			if(a==0.0 && b==0.0) {
				sum->ReleaseParameter(n);
			} else {
				sum->SetParLimits(n, a, b);
			}
			++n;
	}

}


void test() {

	const char*  hist_name = "ggg_IM";
	const double r_min = 700;
	const double r_max = 850;
	const int    npx   = 500;
	const double omega_mass     = 782.0;
	const double expected_width =  15.0;


	TH1* h = NULL;
	gDirectory->GetObject(hist_name, h);

	if(!h) {
		cerr << "Can't find histogram" << endl;
		return;
	}



	TF1* sig = new TF1("sig", "gaus", r_min, r_max);
	sig->SetLineColor(kGreen);
	sig->SetNpx(npx);

	// height
	sig->SetParameter(0, 0.5 * h->GetMaximum());

	// position
	sig->FixParameter(1, omega_mass);

	// width
	sig->SetParameter(2, expected_width);
	


	TF1* bg = new TF1("bg", "pol2", r_min, r_max);
	bg->SetLineColor(kBlue);
	bg->SetNpx(npx);

	bg->SetParameter(0,0);
	bg->SetParName(0, "BG p_{0}");
	bg->SetParameter(1,0);
	bg->SetParName(1, "BG p_{1}");
	bg->SetParameter(2,0);
	bg->SetParName(2, "BG p_{2}");


	TF1* sum = new TF1("sum", "bg+sig", r_min, r_max);
	sum->SetLineColor(kBlack);
	sum->SetNpx(npx);
	SumSetNames(sum, sig, bg);

	CopyParLimits(sum, sig, bg);
	
	sum->Print();

	for(int i=0; i<sum->GetNpar(); ++i){
		double a,b;
		sum->GetParLimits(i,a,b);
		cout << i << " " << a << " " << b << endl;
	}


	TCanvas* c = new TCanvas();
	c->SetTitle(Form("Fit to %s", hist_name));
	h->SetStats(true);
	h->Draw();
	h->Fit(sum, "REM0NB");
	//TODO: Find a way to switch on SetOptFit(). This is a TPaveStats option.

	CopyBack(sum, sig, bg);

	sum->Draw("same");
	bg->Draw("same");


	const double total_area = sum->Integral(r_min, r_max);
	const double bg_area    =  bg->Integral(r_min, r_max);
	const double sig_area   = total_area - bg_area;

	const double sig_to_bg = sig_area / bg_area;
	const double peak_pos  = sig->GetParameter(1);

	cout << "Mass offset = " << peak_pos - omega_mass << " MeV\n";
	cout << "Sig/BG      = " << sig_to_bg << "\n";
	cout << "Sig         = " << sig_area << endl;

	// TODO: choose a position. Positions for TLatex are histogram coordinates.
	TLatex* label = new TLatex(r_min, h->GetMaximum(),Form("Signal content = %lf", sig_area));
	label->Draw();

}

