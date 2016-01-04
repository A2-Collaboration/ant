// use this macro with root as follows:
// root Ant_CBTaggTAPS_5711.dat.xz.root 'rename.C("TAPS_ShortEnergy/RelativeGains","TAPS_ShortEnergy/rel_gamma")'
// mind the shell-escaping of the arguments to rename.C (usually single quotes do the job)

void rename(const char* oldname, const char* newname)
{
	_file0->ReOpen("UPDATE");

	TObject* o = _file0->Get(oldname);
	if(!o) {
		cerr << "Object '" << oldname << "' not found in _file0" << endl;
		exit(1);
	}

	TNamed* n = dynamic_cast<TNamed*>(o);
	if(!n) {
		cerr << "Object '" << oldname << "' does not derive from TNamed" << endl;
		exit(1);
	}

	// this simply adds another TKey (I guess),
	// the object is still there with its old name
	n->SetName(newname);
	_file0->Write();

	exit(0);
}
