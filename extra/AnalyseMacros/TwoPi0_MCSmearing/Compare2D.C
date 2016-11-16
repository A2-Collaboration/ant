// root <mc.root.a.root> <data.root.a.root> Compate2D.C
{
	TDirectory* dmc;
	_file0->GetObject("ETheta", dmc);


	TDirectory* dd;
	_file1->GetObject("ETheta", dd);

	ant::TwoPi0_MCSmearing_Tool::CompareMCData2D(_file0, _file1, "ETheta");

	TH2* sigma;
	gDirectory->GetObject("datamc_simga", sigma);
	ant::TwoPi0_MCSmearing_Tool::getInspectorCanvas(sigma, "proj", dmc ,"MC",dd, "Data");
}
