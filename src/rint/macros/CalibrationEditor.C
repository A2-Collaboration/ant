{
	if ( gSystem->Load("/opt/ant/build/lib/librint.so") != 0 )
	{
		cerr << "  Could not load ant's librint!!!" << endl;
		gSystem->Exit(1);
	}

	ant::TCalibrationEditor e;

	int retval = 0;
	TCanvas* can = new TCanvas("ACE","Ant Calibration Editor");
	while (retval == 0)
	{
		cout << endl;
		e.ListCommands();
	    cout << endl << "-->  ";
		string command;
		cin >> command;
		can->cd();
		retval = e.Execute(command);
		can->Draw();
		can->Update();
	}

	return 0;
}
