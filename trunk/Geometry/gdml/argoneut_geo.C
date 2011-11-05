typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

void argoneut_geo(TString volName="volCryostat"){

gSystem->Load("libGeom");
gSystem->Load("libGdml");

TGeoManager::Import("argoneut.gdml");

drawopt optArgoNeuT[] = {
//   {"volWorld",                 0},
//   {"volDetEnclosure",          kWhite},
//   {"volCryostat",              kOrange},
//   {"volTPCWirePlaneLengthSide", kCyan+3},
//   {"volTPCWirePlaneWidthSide", kRed},
  {"volTPCWire45", kRed},
  {"volTPCWire0", kBlue},
  {"volTPCWidthFace", kRed},
  {"volTPCLengthFace", kCyan+5},
  {"volTPCBottomFace", kOrange},
{"volTubLongSide", kOrange+7},
{"volTubShortSide", kOrange+7},
{"volTubBottom", kOrange+7},
{"voltheX", kOrange+7},
{"volTPCShieldPlane", kBlue},
{"volArgon_solid_L", kRed},
{"volArgon_cap_L", kOrange},
{"volDetEnclosure", kBlue},
{"volMND", kBlue},
{"volTubLongSide", kBlue},
{"volTPCActive",kGreen},
  {0, 0}
};

for (int i=0;; ++i) {
  if (optArgoNeuT[i].volume==0) break;
    gGeoManager->FindVolumeFast(optArgoNeuT[i].volume)->SetLineColor(optArgoNeuT[i].color);
}

TList* mat = gGeoManager->GetListOfMaterials();
TIter next(mat);
TObject *obj;
while (obj = next()) {
 obj->Print();
}

 gGeoManager->CheckOverlaps(0.01);
 gGeoManager->PrintOverlaps();
 gGeoManager->SetMaxVisNodes(70000);

 //gGeoManager->GetTopVolume()->Draw();
 //gGeoManager->FindVolumeFast(volName)->Draw();

 TFile *tf = new TFile("argoneut.root", "RECREATE");
 gGeoManager->Write();
 tf->Close();
}
