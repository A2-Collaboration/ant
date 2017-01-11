//Create a histogram of the CB resolution (MeV) calculated by the official formaula
TH2* CbNomERes() {
	TH2* h = new TH2D("","",16,0,1600,35,-1,1);
	for(int x = 1; x<=16; x++) {
		const double E = h->GetXaxis()->GetBinCenter(x);
		for(y=1;y<=35;y++) {
//			const double t = h->GetYaxis()->GetBinCenter(y);
			const double res = 0.02*E/(TMath::Power(E/1000.0,0.25));
			h->SetBinContent(x,y,res);
		}
	}
	return h;
}

//Divide a histogram by the CB resolution (MeV) calculated by the official formaula
void DivNormCB(TH2* h) {
	for(int x = 1; x<=h->GetNbinsX(); x++) {
		const double E = h->GetXaxis()->GetBinCenter(x);
		for(int y=1;y<=h->GetNbinsY();y++) {
//			const double t = h->GetYaxis()->GetBinCenter(y);
			const double res = 0.02*E/(TMath::Power(E/1000.0,0.25));
			h->SetBinContent(x,y,h->GetBinContent(x,y)/res);
		}
	}
}
