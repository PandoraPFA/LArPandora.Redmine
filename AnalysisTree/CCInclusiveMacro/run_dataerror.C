{
gROOT->Reset();
gROOT->LoadMacro("dataerror.C");
// gROOT->LoadMacro("bgkin.C");
datakin t;
t->Loop();
}
