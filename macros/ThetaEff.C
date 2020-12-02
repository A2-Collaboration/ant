void PlotStuff()
{
	UInt_t n;
	TString name;

	TCanvas *c1 = new TCanvas( "canvas", "Testing", 800, 0, 700, 700);
	c1->Divide( 2, 2);

	name = "../out/He4Pi0_300MeV.root";
	TFile *geant = new TFile( name);

	name = "He4Pi0/h3D_MEpi0";
	TH3D *h3 = (TH3D*) geant->Get( name);
	h3->GetXaxis()->SetRangeUser( -100, 10);
//	h3->GetXaxis()->SetRangeUser( -20, -10);

	TH1D *histx = (TH1D*) h3->Project3D( "x");
	TH1D *histy = (TH1D*) h3->Project3D( "y");
	TH1D *hyclone = (TH1D*)histy->Clone( "hyclone");

//	n = histy->GetXaxis()->GetNbins();
//	cout << n << endl;

	c1->cd( 1);
	histx->Draw();
	c1->cd( 2);
	hyclone->Draw();

	name = "../evgen/5cm/pi0_he4_300_in.root";
	TFile *evgen = new TFile( name);
	name = "h4";
	TH1D *ev_hist = (TH1D*) evgen->Get( name);
	c1->cd( 3);
	ev_hist->Rebin( 10);
	ev_hist->Draw();

	histy->Sumw2();
	histy->Divide( ev_hist);
	c1->cd( 4);
	histy->Draw();

//	n = ev_hist->GetXaxis()->GetNbins();
//	cout << n << endl;

}

void ThetaEff( UInt_t rebin = 1, Bool_t plot = kTRUE)
{
	UInt_t i, end;
	Double_t yld_ge, yld_ev, eff, deff;
	TString name;

	name = "compton_p_200_in.root";
	TFile *evgen = new TFile( name);

	name = "ARH_Compton.root";
	TFile *geant = new TFile( name);

	name = "h4";
	TH1D *ev_hist = (TH1D*) evgen->Get( name);
	TH1D *ev_clone = (TH1D*)ev_hist->Clone( "evgen");

	name = "PHYS_PhotonTheta";
	TH1D *h1 = (TH1D*) geant->Get( name);
	TH1D *geant_clone = (TH1D*)h1->Clone( "geant");

	yld_ge = h1->Integral();
	yld_ev = ev_hist->Integral();
	eff = (double) yld_ge/yld_ev;
	deff = sqrt( yld_ge/(yld_ev*yld_ev) + (yld_ge*yld_ge)/pow( yld_ev, 3));

	ev_hist->Rebin( rebin);
	h1->Rebin( rebin);
	h1->Sumw2();
	h1->Divide( ev_hist);

	name = Form( "  0 %6.4f %6.4f", eff, deff);
	cout << name << endl;

	end = 180/rebin;

	for ( i = 1; i <= end; i++) {
		name = Form( "%3.0f %6.4f %6.4f", h1->GetBinCenter( i),
				h1->GetBinContent( i), h1->GetBinError( i));
		cout << name << endl;
	}

	if ( plot == kTRUE) {
		TCanvas *c1 = new TCanvas( "canvas", "Set-Up Efficiency", 700, 700);
		c1->Divide(1,2);

		c1->cd(1);
		ev_clone->Draw();
		ev_clone->SetFillColor( 46);
		ev_clone->SetLineWidth( 2);
		ev_clone->SetMinimum( 0);
		geant_clone->Draw( "same");
		geant_clone->SetFillColor( 38);
		geant_clone->SetLineWidth( 2);

		c1->cd(2);
		h1->Draw( "E");
		h1->SetLineWidth( 2);
		h1->SetMinimum( 0);
		h1->SetMaximum( 1.5);
	}
}

void Testing()
{
	TString name;

	name = "ARH_Compton.root";
	TFile *geant = new TFile( name);

	name = "PHYS_MissingMassPrompt_v_PhotonThetaPrompt";
	TH2D *hist2D = (TH2D*) geant->Get( name);

	hist2D->Draw();

}

void Testing2( UInt_t nphotbin = 0, UInt_t rebin = 1)
{
	TString name;

	TCanvas *c1 = new TCanvas( "canvas", "Set-Up Efficiency", 1250, 0, 800, 1200);
	c1->Divide( 1, 2);

	name = "ARH_Compton.root";
	TFile *geant = new TFile( name);

	name = "PHYS_NPhoton_v_PhotonThetaPrompt";
	TH2D *hist2D = (TH2D*) geant->Get( name);

	if ( nphotbin != 0) hist2D->GetYaxis()->SetRange( nphotbin, nphotbin);
	TH1D *histx = hist2D->ProjectionX( "theta_eff_cut");

	name = "compton_p_200_in.root";
	TFile *evgen = new TFile( name);
	name = "h4";
	TH1D *ev_hist = (TH1D*) evgen->Get( name);
	TH1D *ev_clone = (TH1D*)ev_hist->Clone( "evgen");

	TH1D *geant_clone = (TH1D*)histx->Clone( "geant");

	ev_clone->Rebin( rebin);
	geant_clone->Rebin( rebin);
	geant_clone->Sumw2();
	geant_clone->Divide( ev_clone);

	c1->cd( 2);
	geant_clone->SetMaximum( 1.2);
	geant_clone->Draw();

	c1->cd( 1);
	ev_hist->Draw();
//	ev_hist->SetFillColor( 46);
	ev_hist->SetLineWidth( 4);
	histx->Draw( "same");
//	histx->SetFillColor( 38);
//	histx->SetLineStyle( 2);
	histx->SetLineWidth( 4);
	histx->SetLineColor( 2);

}
