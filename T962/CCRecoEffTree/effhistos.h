//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 20 12:38:05 2012 by ROOT version 5.30/02
// from TTree EffTree/EffTree
// found on file: ccrecoefftree.root
//////////////////////////////////////////////////////////

#ifndef effhistos_h
#define effhistos_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class effhistos {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Double_t        pot;
   Int_t           isdata;
   Double_t        vtxx_reco;
   Double_t        vtxy_reco;
   Double_t        vtxz_reco;
   Double_t        trackstart_dcosx_reco;
   Double_t        trackstart_dcosy_reco;
   Double_t        trackstart_dcosz_reco;
   Double_t        trackexit_dcosx_reco;
   Double_t        trackexit_dcosy_reco;
   Double_t        trackexit_dcosz_reco;
   Double_t        trackstart_x_reco;
   Double_t        trackstart_y_reco;
   Double_t        trackstart_z_reco;
   Double_t        trackexit_x_reco;
   Double_t        trackexit_y_reco;
   Double_t        trackexit_z_reco;
   Int_t           nmatched_reco;
   Double_t        trk_mom_minos;
   Double_t        trk_charge_minos;
   Double_t        trk_dcosx_minos;
   Double_t        trk_dcosy_minos;
   Double_t        trk_dcosz_minos;
   Double_t        trk_vtxx_minos;
   Double_t        trk_vtxy_minos;
   Double_t        trk_vtxz_minos;
   Int_t           test_charge_minos;
   Int_t           nuPDG_truth;
   Int_t           ccnc_truth;
   Int_t           mode_truth;
   Double_t        enu_truth;
   Double_t        Q2_truth;
   Int_t           hitnuc_truth;
   Double_t        W_truth;
   Double_t        nuvtxx_truth;
   Double_t        nuvtxy_truth;
   Double_t        nuvtxz_truth;
   Double_t        lep_mom_truth;
   Double_t        lep_dcosx_truth;
   Double_t        lep_dcosy_truth;
   Double_t        lep_dcosz_truth;
   Int_t           nminos_tracks;
   Double_t        trk_charge_minos_all[1];   //[nminos_tracks]
   Int_t           muon_reco;
   Int_t           minos_enter_true;
   Double_t        trackstart_x_reco_muon;
   Double_t        trackstart_y_reco_muon;
   Double_t        trackstart_z_reco_muon;
   Double_t        trackexit_x_reco_muon;
   Double_t        trackexit_y_reco_muon;
   Double_t        trackexit_z_reco_muon;
   Int_t           muon_exits;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_pot;   //!
   TBranch        *b_isdata;   //!
   TBranch        *b_vtxx_reco;   //!
   TBranch        *b_vtxy_reco;   //!
   TBranch        *b_vtxz_reco;   //!
   TBranch        *b_trackstart_dcosx_reco;   //!
   TBranch        *b_trackstart_dcosy_reco;   //!
   TBranch        *b_trackstart_dcosz_reco;   //!
   TBranch        *b_trackexit_dcosx_reco;   //!
   TBranch        *b_trackexit_dcosy_reco;   //!
   TBranch        *b_trackexit_dcosz_reco;   //!
   TBranch        *b_trackstart_x_reco;   //!
   TBranch        *b_trackstart_y_reco;   //!
   TBranch        *b_trackstart_z_reco;   //!
   TBranch        *b_trackexit_x_reco;   //!
   TBranch        *b_trackexit_y_reco;   //!
   TBranch        *b_trackexit_z_reco;   //!
   TBranch        *b_nmatched_reco;   //!
   TBranch        *b_trk_mom_minos;   //!
   TBranch        *b_trk_charge_minos;   //!
   TBranch        *b_trk_dcosx_minos;   //!
   TBranch        *b_trk_dcosy_minos;   //!
   TBranch        *b_trk_dcosz_minos;   //!
   TBranch        *b_trk_vtxx_minos;   //!
   TBranch        *b_trk_vtxy_minos;   //!
   TBranch        *b_trk_vtxz_minos;   //!
   TBranch        *b_test_charge_minos;   //!
   TBranch        *b_nuPDG_truth;   //!
   TBranch        *b_ccnc_truth;   //!
   TBranch        *b_mode_truth;   //!
   TBranch        *b_enu_truth;   //!
   TBranch        *b_Q2_truth;   //!
   TBranch        *b_hitnuc_truth;   //!
   TBranch        *b_W_truth;   //!
   TBranch        *b_nuvtxx_truth;   //!
   TBranch        *b_nuvtxy_truth;   //!
   TBranch        *b_nuvtxz_truth;   //!
   TBranch        *b_lep_mom_truth;   //!
   TBranch        *b_lep_dcosx_truth;   //!
   TBranch        *b_lep_dcosy_truth;   //!
   TBranch        *b_lep_dcosz_truth;   //!
   TBranch        *b_nminos_tracks;   //!
   TBranch        *b_trk_charge_minos_all;   //!
   TBranch        *b_muon_reco;   //!
   TBranch        *b_minos_enter_true;   //!
   TBranch        *b_trackstart_x_reco_muon;   //!
   TBranch        *b_trackstart_y_reco_muon;   //!
   TBranch        *b_trackstart_z_reco_muon;   //!
   TBranch        *b_trackexit_x_reco_muon;   //!
   TBranch        *b_trackexit_y_reco_muon;   //!
   TBranch        *b_trackexit_z_reco_muon;   //!
   TBranch        *b_muon_exits;   //!

   effhistos(TTree *tree=0);
   virtual ~effhistos();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef effhistos_cxx
effhistos::effhistos(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ccrecoefftree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ccrecoefftree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ccrecoefftree.root:/ccrecoefftree");
      dir->GetObject("EffTree",tree);

   }
   Init(tree);
}

effhistos::~effhistos()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t effhistos::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t effhistos::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void effhistos::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("pot", &pot, &b_pot);
   fChain->SetBranchAddress("isdata", &isdata, &b_isdata);
   fChain->SetBranchAddress("vtxx_reco", &vtxx_reco, &b_vtxx_reco);
   fChain->SetBranchAddress("vtxy_reco", &vtxy_reco, &b_vtxy_reco);
   fChain->SetBranchAddress("vtxz_reco", &vtxz_reco, &b_vtxz_reco);
   fChain->SetBranchAddress("trackstart_dcosx_reco", &trackstart_dcosx_reco, &b_trackstart_dcosx_reco);
   fChain->SetBranchAddress("trackstart_dcosy_reco", &trackstart_dcosy_reco, &b_trackstart_dcosy_reco);
   fChain->SetBranchAddress("trackstart_dcosz_reco", &trackstart_dcosz_reco, &b_trackstart_dcosz_reco);
   fChain->SetBranchAddress("trackexit_dcosx_reco", &trackexit_dcosx_reco, &b_trackexit_dcosx_reco);
   fChain->SetBranchAddress("trackexit_dcosy_reco", &trackexit_dcosy_reco, &b_trackexit_dcosy_reco);
   fChain->SetBranchAddress("trackexit_dcosz_reco", &trackexit_dcosz_reco, &b_trackexit_dcosz_reco);
   fChain->SetBranchAddress("trackstart_x_reco", &trackstart_x_reco, &b_trackstart_x_reco);
   fChain->SetBranchAddress("trackstart_y_reco", &trackstart_y_reco, &b_trackstart_y_reco);
   fChain->SetBranchAddress("trackstart_z_reco", &trackstart_z_reco, &b_trackstart_z_reco);
   fChain->SetBranchAddress("trackexit_x_reco", &trackexit_x_reco, &b_trackexit_x_reco);
   fChain->SetBranchAddress("trackexit_y_reco", &trackexit_y_reco, &b_trackexit_y_reco);
   fChain->SetBranchAddress("trackexit_z_reco", &trackexit_z_reco, &b_trackexit_z_reco);
   fChain->SetBranchAddress("nmatched_reco", &nmatched_reco, &b_nmatched_reco);
   fChain->SetBranchAddress("trk_mom_minos", &trk_mom_minos, &b_trk_mom_minos);
   fChain->SetBranchAddress("trk_charge_minos", &trk_charge_minos, &b_trk_charge_minos);
   fChain->SetBranchAddress("trk_dcosx_minos", &trk_dcosx_minos, &b_trk_dcosx_minos);
   fChain->SetBranchAddress("trk_dcosy_minos", &trk_dcosy_minos, &b_trk_dcosy_minos);
   fChain->SetBranchAddress("trk_dcosz_minos", &trk_dcosz_minos, &b_trk_dcosz_minos);
   fChain->SetBranchAddress("trk_vtxx_minos", &trk_vtxx_minos, &b_trk_vtxx_minos);
   fChain->SetBranchAddress("trk_vtxy_minos", &trk_vtxy_minos, &b_trk_vtxy_minos);
   fChain->SetBranchAddress("trk_vtxz_minos", &trk_vtxz_minos, &b_trk_vtxz_minos);
   fChain->SetBranchAddress("test_charge_minos", &test_charge_minos, &b_test_charge_minos);
   fChain->SetBranchAddress("nuPDG_truth", &nuPDG_truth, &b_nuPDG_truth);
   fChain->SetBranchAddress("ccnc_truth", &ccnc_truth, &b_ccnc_truth);
   fChain->SetBranchAddress("mode_truth", &mode_truth, &b_mode_truth);
   fChain->SetBranchAddress("enu_truth", &enu_truth, &b_enu_truth);
   fChain->SetBranchAddress("Q2_truth", &Q2_truth, &b_Q2_truth);
   fChain->SetBranchAddress("hitnuc_truth", &hitnuc_truth, &b_hitnuc_truth);
   fChain->SetBranchAddress("W_truth", &W_truth, &b_W_truth);
   fChain->SetBranchAddress("nuvtxx_truth", &nuvtxx_truth, &b_nuvtxx_truth);
   fChain->SetBranchAddress("nuvtxy_truth", &nuvtxy_truth, &b_nuvtxy_truth);
   fChain->SetBranchAddress("nuvtxz_truth", &nuvtxz_truth, &b_nuvtxz_truth);
   fChain->SetBranchAddress("lep_mom_truth", &lep_mom_truth, &b_lep_mom_truth);
   fChain->SetBranchAddress("lep_dcosx_truth", &lep_dcosx_truth, &b_lep_dcosx_truth);
   fChain->SetBranchAddress("lep_dcosy_truth", &lep_dcosy_truth, &b_lep_dcosy_truth);
   fChain->SetBranchAddress("lep_dcosz_truth", &lep_dcosz_truth, &b_lep_dcosz_truth);
   fChain->SetBranchAddress("nminos_tracks", &nminos_tracks, &b_nminos_tracks);
   fChain->SetBranchAddress("trk_charge_minos_all", &trk_charge_minos_all, &b_trk_charge_minos_all);
   fChain->SetBranchAddress("muon_reco", &muon_reco, &b_muon_reco);
   fChain->SetBranchAddress("minos_enter_true", &minos_enter_true, &b_minos_enter_true);
   fChain->SetBranchAddress("trackstart_x_reco_muon", &trackstart_x_reco_muon, &b_trackstart_x_reco_muon);
   fChain->SetBranchAddress("trackstart_y_reco_muon", &trackstart_y_reco_muon, &b_trackstart_y_reco_muon);
   fChain->SetBranchAddress("trackstart_z_reco_muon", &trackstart_z_reco_muon, &b_trackstart_z_reco_muon);
   fChain->SetBranchAddress("trackexit_x_reco_muon", &trackexit_x_reco_muon, &b_trackexit_x_reco_muon);
   fChain->SetBranchAddress("trackexit_y_reco_muon", &trackexit_y_reco_muon, &b_trackexit_y_reco_muon);
   fChain->SetBranchAddress("trackexit_z_reco_muon", &trackexit_z_reco_muon, &b_trackexit_z_reco_muon);
   fChain->SetBranchAddress("muon_exits", &muon_exits, &b_muon_exits);
   Notify();
}

Bool_t effhistos::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void effhistos::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t effhistos::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef effhistos_cxx
