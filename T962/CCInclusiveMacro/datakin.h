//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  5 19:16:41 2011 by ROOT version 5.28/00c
// from TTree anatree/analysis tree
// found on file: simkinreco_hist_iso.root
//////////////////////////////////////////////////////////

#ifndef datakin_h
#define datakin_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class datakin {
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
   Int_t           nclusu_reco;
   Int_t           nclusv_reco;
   Int_t           nclusw_reco;
   Int_t           ntracks_reco;
   Int_t           nvertextracks_reco;
   Int_t           nvertexclustersu_reco;
   Int_t           nvertexclustersv_reco;
   Int_t           nvertexclustersw_reco;
   Int_t           ntrackendonboundary_reco;
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
   Float_t         mc_index_minos;
   Double_t        mc_pdg_minos;
   Double_t        mc_px_minos;
   Double_t        mc_py_minos;
   Double_t        mc_pz_minos;
   Double_t        mc_ene_minos;
   Double_t        mc_mass_minos;
   Double_t        mc_vtxx_minos;
   Double_t        mc_vtxy_minos;
   Double_t        mc_vtxz_minos;
   Int_t           trkcontained_minos;
   Int_t           test_charge_minos;
   Double_t        vtxx_scan;
   Double_t        vtxy_scan;
   Double_t        vtxz_scan;
   Int_t           ntracks_scan;
   Int_t           nshowers_scan;
   Int_t           neutrino_scan;
   Int_t           maybeneutrino_scan;
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

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_pot;   //!
   TBranch        *b_isdata;   //!
   TBranch        *b_vtxx_reco;   //!
   TBranch        *b_vtxy_reco;   //!
   TBranch        *b_vtxz_reco;   //!
   TBranch        *b_nclusu_reco;   //!
   TBranch        *b_nclusv_reco;   //!
   TBranch        *b_nclusw_reco;   //!
   TBranch        *b_ntracks_reco;   //!
   TBranch        *b_nvertextracks_reco;   //!
   TBranch        *b_nvertexclustersu_reco;   //!
   TBranch        *b_nvertexclustersv_reco;   //!
   TBranch        *b_nvertexclustersw_reco;   //!
   TBranch        *b_ntrackendonboundary_reco;   //!
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
   TBranch        *b_mc_index_minos;   //!
   TBranch        *b_mc_pdg_minos;   //!
   TBranch        *b_mc_px_minos;   //!
   TBranch        *b_mc_py_minos;   //!
   TBranch        *b_mc_pz_minos;   //!
   TBranch        *b_mc_ene_minos;   //!
   TBranch        *b_mc_mass_minos;   //!
   TBranch        *b_mc_vtxx_minos;   //!
   TBranch        *b_mc_vtxy_minos;   //!
   TBranch        *b_mc_vtxz_minos;   //!
   TBranch        *b_trkcontained_minos;   //!
   TBranch        *b_test_charge_minos;   //!
   TBranch        *b_vtxx_scan;   //!
   TBranch        *b_vtxy_scan;   //!
   TBranch        *b_vtxz_scan;   //!
   TBranch        *b_ntracks_scan;   //!
   TBranch        *b_nshowers_scan;   //!
   TBranch        *b_neutrino_scan;   //!
   TBranch        *b_maybeneutrino_scan;   //!
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

   datakin(TTree *tree=0);
   virtual ~datakin();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef datakin_cxx
datakin::datakin(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("simkinreco_hist_iso.root");
      if (!f) {
         f = new TFile("simkinreco_hist_iso.root");
         f->cd("simkinreco_hist_iso.root:/analysistree");
      }
      tree = (TTree*)gDirectory->Get("anatree");

   }
   Init(tree);
}

datakin::~datakin()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t datakin::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t datakin::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void datakin::Init(TTree *tree)
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
   fChain->SetBranchAddress("nclusu_reco", &nclusu_reco, &b_nclusu_reco);
   fChain->SetBranchAddress("nclusv_reco", &nclusv_reco, &b_nclusv_reco);
   fChain->SetBranchAddress("nclusw_reco", &nclusw_reco, &b_nclusw_reco);
   fChain->SetBranchAddress("ntracks_reco", &ntracks_reco, &b_ntracks_reco);
   fChain->SetBranchAddress("nvertextracks_reco", &nvertextracks_reco, &b_nvertextracks_reco);
   fChain->SetBranchAddress("nvertexclustersu_reco", &nvertexclustersu_reco, &b_nvertexclustersu_reco);
   fChain->SetBranchAddress("nvertexclustersv_reco", &nvertexclustersv_reco, &b_nvertexclustersv_reco);
   fChain->SetBranchAddress("nvertexclustersw_reco", &nvertexclustersw_reco, &b_nvertexclustersw_reco);
   fChain->SetBranchAddress("ntrackendonboundary_reco", &ntrackendonboundary_reco, &b_ntrackendonboundary_reco);
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
   fChain->SetBranchAddress("mc_index_minos", &mc_index_minos, &b_mc_index_minos);
   fChain->SetBranchAddress("mc_pdg_minos", &mc_pdg_minos, &b_mc_pdg_minos);
   fChain->SetBranchAddress("mc_px_minos", &mc_px_minos, &b_mc_px_minos);
   fChain->SetBranchAddress("mc_py_minos", &mc_py_minos, &b_mc_py_minos);
   fChain->SetBranchAddress("mc_pz_minos", &mc_pz_minos, &b_mc_pz_minos);
   fChain->SetBranchAddress("mc_ene_minos", &mc_ene_minos, &b_mc_ene_minos);
   fChain->SetBranchAddress("mc_mass_minos", &mc_mass_minos, &b_mc_mass_minos);
   fChain->SetBranchAddress("mc_vtxx_minos", &mc_vtxx_minos, &b_mc_vtxx_minos);
   fChain->SetBranchAddress("mc_vtxy_minos", &mc_vtxy_minos, &b_mc_vtxy_minos);
   fChain->SetBranchAddress("mc_vtxz_minos", &mc_vtxz_minos, &b_mc_vtxz_minos);
   fChain->SetBranchAddress("trkcontained_minos", &trkcontained_minos, &b_trkcontained_minos);
   fChain->SetBranchAddress("test_charge_minos", &test_charge_minos, &b_test_charge_minos);
   fChain->SetBranchAddress("vtxx_scan", &vtxx_scan, &b_vtxx_scan);
   fChain->SetBranchAddress("vtxy_scan", &vtxy_scan, &b_vtxy_scan);
   fChain->SetBranchAddress("vtxz_scan", &vtxz_scan, &b_vtxz_scan);
   fChain->SetBranchAddress("ntracks_scan", &ntracks_scan, &b_ntracks_scan);
   fChain->SetBranchAddress("nshowers_scan", &nshowers_scan, &b_nshowers_scan);
   fChain->SetBranchAddress("neutrino_scan", &neutrino_scan, &b_neutrino_scan);
   fChain->SetBranchAddress("maybeneutrino_scan", &maybeneutrino_scan, &b_maybeneutrino_scan);
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
   Notify();
}

Bool_t datakin::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void datakin::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t datakin::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef datakin_cxx
