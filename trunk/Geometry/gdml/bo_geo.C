typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

bo_geo(TString volName="volCryostat"){

  gSystem->Load("libGeom");
  gSystem->Load("libGdml");
  
  TGeoManager::Import("bo.gdml");
  
  drawopt optbo[] = {
    {"volWorld",                 0},
    {"volCryostat",      kRed},
    {"volBoOuterShell",  kOrange},
    {"volBoOuterShellBottom", kOrange},
    {"volBoInnerShell",  kOrange-3},
    {"volBoInnerShellBottom", kOrange-3},
    {"volBoTopLid",      kYellow},
    {"volTPCSheet",        kCyan+3},
    {"volTPCHalfStrip",	 kCyan+1},
    {"volTPCStrip",		kCyan+1},
    {"volTPCSR",			kCyan+2},
    {"volTPCBottomRing",		kCyan-5},
    //  {"volTube",			kOrange},
    //  {"volTubeCorner",		kOrange},
    {"volTPCLROuterPart",	kCyan+2},
    {"volTPCLREH",		kCyan+2},
    {"volTPCLRBoxBetweenEHs",	kCyan+2},
    {"volTPCLRBoxBesideEHs",	kCyan+2},
    {"volTPCLRBoxTangentToEHs",	kCyan+2},
    {"volTPCLRInnerPart",	kCyan+2},
    {"volTPCWire",		kGray},
    {"volTPCRing",		kCyan-5},
    {"volTPCPlane",		kCyan},
    {"volDetEnclosure",		kWhite},
    {0, 0}
  };
  
  for (int i=0;; ++i) {
    if (optbo[i].volume==0) break;
    gGeoManager->FindVolumeFast(optbo[i].volume)->SetLineColor(optbo[i].color);
  }
  
  gGeoManager->CheckOverlaps(0.01);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);
  
  
  TList* mat = gGeoManager->GetListOfMaterials();
  TIter next(mat);
  TObject *obj;
  while (obj = next()) {
    obj->Print();
  }

//gGeoManager->GetTopVolume()->Draw();
  gGeoManager->FindVolumeFast(volName)->Draw();

  TFile *tf = new TFile("bo.root", "RECREATE");
  gGeoManager->Write();
  tf->Close();
}
