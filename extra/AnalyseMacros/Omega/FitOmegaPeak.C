#include "TF1.h"
#include "TH1.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include "root-addons/analysis_codes/Math.h"

using namespace std;
using namespace ant;

void FitOmegaPeak(const bool fixOmegaMass=false, const double r_min=700.0, const double r_max=850.0) {

	const char*  hist_name = "ggg_IM";
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
	if(fixOmegaMass)
		sig->FixParameter(1, omega_mass);
	else
		sig->SetParameter(1, omega_mass);

	// width
	sig->SetParameter(2, expected_width);



	TF1* bg = new TF1("bg", "pol2", r_min, r_max);
	bg->SetLineColor(kBlue);

	bg->SetParameter(0,0);
	bg->SetParName(0, "BG p_{0}");
	bg->SetParameter(1,0);
	bg->SetParName(1, "BG p_{1}");
	bg->SetParameter(2,0);
	bg->SetParName(2, "BG p_{2}");

	TFSum* sum = new TFSum("sum", sig, bg, r_min, r_max);
	sum->SetNpx(npx);


	TCanvas* c = new TCanvas();
	c->SetTitle(Form("Fit to %s", hist_name));
	h->SetStats(true);
	gStyle->SetOptFit(1);
	h->Draw();
	h->Fit(sum->Function(), "REM0NB");

	sum->SyncToFcts();

	sum->Draw();


	const double total_area = sum->Function()->Integral(r_min, r_max);
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

