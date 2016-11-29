// root <mc.root.a.root> <data.root.a.root> Compate2D.C
void Compare2D(const char* folder) {
	TDirectory* dmc;
	_file0->GetObject(folder, dmc);


	TDirectory* dd;
	_file1->GetObject(folder, dd);

	ant::TwoPi0_MCSmearing_Tool::CompareMCData2D(_file0, _file1, folder);

	TH2* sigma;
	gDirectory->GetObject("datamc_simga", sigma);
	ant::TwoPi0_MCSmearing_Tool::getInspectorCanvas(sigma, "proj", dmc ,"MC",dd, "Data");
}
