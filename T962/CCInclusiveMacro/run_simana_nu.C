{
	gROOT->Reset();
	gROOT->LoadMacro("simana_nu.C");
	simkinana t;
	t->Loop();
}
