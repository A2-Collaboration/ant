{
	cout << "Ants librint.so loaded correctly: " << ( gSystem->Load("/opt/ant/build/lib/librint.so") == 0 ) << endl << endl;

	ant::TCalibrationEditor e;
	e.AddSomeRandomData();

	e.ShowHistory("testID");
}

