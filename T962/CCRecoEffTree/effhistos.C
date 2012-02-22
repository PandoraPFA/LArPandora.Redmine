#define effhistos_cxx
#include "effhistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void effhistos::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L effhistos.C
//      Root > effhistos t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  double RDIFF_CUT = 27.;//cm
  double THETADIFF_CUT = 0.4;//rad
  double FIDX_CUT = 3.;//cm
  double FIDY_CUT = 4.;//cm
  double FIDZup_CUT = 6.;//cm
  double FIDZdown_CUT = 4.;//cm  
  
  gStyle->SetOptStat(kFALSE);
  
  /// muon momentum histograms 
  TH1D* histeff1111D = new TH1D("", "",20,0.,25);
  histeff1111D->Sumw2();
  
  TH1D* histeff2222 = new TH1D("", "",20,0.,25);
  histeff2222->Sumw2();


  TH1D* histeff1D = new TH1D("", "",20,0.,25);
  histeff1D->Sumw2();
  TH1D* histeff11 = new TH1D("", "",20,0.,25);
  histeff11->Sumw2();

  TH1D* histeff2D = new TH1D("", "",20,0.,25);
  histeff2D->Sumw2();
  TH1D* histeff22 = new TH1D("", "",20,0.,25);
  histeff22->Sumw2();
  
  TH1D* histeff3D = new TH1D("", "",20,0.,25);
  histeff3D->Sumw2();
  TH1D* histeff33 = new TH1D("", "",20,0.,25);
  histeff33->Sumw2();

  TH1D* histeff4D = new TH1D("", "",20,0.,25);
  histeff4D->Sumw2();
  TH1D* histeff44 = new TH1D("", "",20,0.,25);
  histeff44->Sumw2();

  TH1D* histeff55 = new TH1D("", "",20,0.,25);
  histeff55->Sumw2();

  TH1D* histeff66 = new TH1D("", "",20,0.,25);
  histeff66->Sumw2();

 TH1D* histeff0D = new TH1D("", "",20,0.,25);
  histeff0D->Sumw2();

 TH1D* histeff00 = new TH1D("", "",20,0.,25);
  histeff00->Sumw2();

  TH1D* histeffXX = new TH1D("", "",20,0.,25);
  histeffXX->Sumw2();

  TH1D* histeffXD = new TH1D("", "",20,0.,25);
  histeffXD->Sumw2();


/////////////////////////////////////
  /// muon angle histograms
  TH1D* histeffT1111D = new TH1D("", "",18,0.,36);
  histeffT1111D->Sumw2();

  TH1D* histeffT2222 = new TH1D("", "",18,0.,36);
  histeffT2222->Sumw2();

  TH1D* histeffT1D = new TH1D("", "",18,0.,36);
  histeffT1D->Sumw2();
  TH1D* histeffT11 = new TH1D("", "",18,0.,36);
  histeffT11->Sumw2();

  TH1D* histeffT2D = new TH1D("", "",18,0.,36);
  histeffT2D->Sumw2();
  TH1D* histeffT22 = new TH1D("", "",18,0.,36);
  histeffT22->Sumw2();

  TH1D* histeffT3D = new TH1D("", "",18,0.,36);
  histeffT3D->Sumw2();
  TH1D* histeffT33 = new TH1D("", "",18,0.,36);
  histeffT33->Sumw2();
 
  TH1D* histeffT4D = new TH1D("", "",18,0.,36);
  histeffT4D->Sumw2();
  TH1D* histeffT44 = new TH1D("", "",18,0.,36);
  histeffT44->Sumw2();

  TH1D* histeffT55 = new TH1D("", "",18,0.,36);
  histeffT55->Sumw2();

  TH1D* histeffT66 = new TH1D("", "",18,0.,36);
  histeffT66->Sumw2();

  TH1D* histeffT0D = new TH1D("", "",18,0.,36);
  histeffT0D->Sumw2();
  
  TH1D* histeffT00 = new TH1D("", "",18,0.,36);
  histeffT00->Sumw2();

  TH1D* histeffTXX = new TH1D("", "",18,0.,36);
  histeffTXX->Sumw2();

  TH1D* histeffTXD = new TH1D("", "",18,0.,36);
  histeffTXD->Sumw2();


////////////////////////////////////////
  int a = 0; int b = 0; int c = 0; int d = 0; int e = 0; int f = 0; 
  int g = 0; int j = 0; int k = 0; int true = 0; 

  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double dcosztruth=(lep_dcosz_truth*TMath::Cos(.0583497))-(lep_dcosy_truth*TMath::Sin(.0583497));

      if(ccnc_truth==0 && nuPDG_truth==14){
	if(	 
	   nuvtxx_truth>0.+FIDX_CUT 
	   && nuvtxx_truth<47.-FIDX_CUT 
	   && nuvtxy_truth>-20.+FIDY_CUT 
	   && nuvtxy_truth<20.-FIDY_CUT 
	   && nuvtxz_truth>0.+FIDZup_CUT 
	   && nuvtxz_truth<90.-FIDZdown_CUT){
	  
	  histeff0D->Fill(lep_mom_truth);
	  histeffT0D->Fill((180./3.14159)*TMath::ACos(dcosztruth));
	  if(lep_mom_truth>1)
	    true++; // number of true events in FV
	}
	
	if(
	   trackstart_x_reco_muon>0.+FIDX_CUT 
	   &&  trackstart_x_reco_muon<47.-FIDX_CUT  
	   &&  trackstart_y_reco_muon>-20.+FIDY_CUT 
	   &&  trackstart_y_reco_muon<20.-FIDY_CUT
	   &&  trackstart_z_reco_muon>0.+FIDZup_CUT 
	   &&  trackstart_z_reco_muon<90.-FIDZdown_CUT 
	   && muon_reco == 1){ 

	  histeff1D->Fill(lep_mom_truth);
	  histeff00->Fill(lep_mom_truth);
	  
	  histeffT1D->Fill((180./3.14159)*TMath::ACos(dcosztruth));
	  histeffT00->Fill((180./3.14159)*TMath::ACos(dcosztruth));
	  if(lep_mom_truth>1)
	    a++;  // number of cc nu-mu events with a reco muon starting in FV
	  
	  
	  if( vtxx_reco>0.+FIDX_CUT 
	      &&  vtxx_reco<47.-FIDX_CUT  
	      &&  vtxy_reco>-20.+FIDY_CUT 
	      &&  vtxy_reco<20.-FIDY_CUT 
	      &&  vtxz_reco>0.+FIDZup_CUT
	      &&  vtxz_reco<90.-FIDZdown_CUT){
	   
	    histeff11->Fill(lep_mom_truth);
	    histeff2D->Fill(lep_mom_truth);
	    
	    histeffT11->Fill((180./3.14159)*TMath::ACos(dcosztruth));
	    histeffT2D->Fill((180./3.14159)*TMath::ACos(dcosztruth));
	    if(lep_mom_truth>1)
	      b++; // number of cc nu-mu events with a reco muon starting in FV + reco vertex in FV 
	    
	    if(muon_exits == 1){
	      histeff22->Fill(lep_mom_truth);
	      histeff3D->Fill(lep_mom_truth);
	      
	      histeffT22->Fill((180./3.14159)*TMath::ACos(dcosztruth));
	      histeffT3D->Fill((180./3.14159)*TMath::ACos(dcosztruth));  
	      if(lep_mom_truth>1)
		c++; // number of cc nu-mu events with a reco muon starting in FV + reco vertex in FV + muon exits the TPC  
	    	      
	      if(nminos_tracks>0){
		
		histeff33->Fill(lep_mom_truth);
		histeff4D->Fill(lep_mom_truth);
		
		histeffT33->Fill((180./3.14159)*TMath::ACos(dcosztruth));
		histeffT4D->Fill((180./3.14159)*TMath::ACos(dcosztruth));
		if(lep_mom_truth>1)
		  d++; // number of cc nu-mu events with a reco muon starting in FV + reco vertex in FV + muon exits the TPC + minos tracks>0
		
		if(test_charge_minos==-1){ 
		  histeffXX->Fill(lep_mom_truth);
		  histeffXD->Fill(lep_mom_truth);

		  histeffTXX->Fill((180./3.14159)*TMath::ACos(dcosztruth));
		  histeffTXD->Fill((180./3.14159)*TMath::ACos(dcosztruth));		  
		  if(lep_mom_truth>1)
		    e++; // number of cc nu-mu events with a reco muon starting in FV + reco vertex in FV + muon exits the TPC + minos tracks>0 + atleast 1 -vely charged MINOS track

		  if(nmatched_reco == 1 && trk_charge_minos<0 && trk_mom_minos>0){
		    histeff44->Fill(lep_mom_truth);
		    
		    histeffT44->Fill((180./3.14159)*TMath::ACos(dcosztruth));
		    if(lep_mom_truth>1)
		      f++; // number of cc nu-mu events with a reco muon starting in FV + reco vertex in FV + muon exits the TPC + minos tracks>0 + atleast 1 -vely charged MINOS track + matched with correct sign
		    
		    double lardirectionEnd1[3]={ trackexit_dcosx_reco,  trackexit_dcosy_reco,  trackexit_dcosz_reco};
		    double larEnd1[3]= {  trackexit_x_reco,  trackexit_y_reco,  trackexit_z_reco};
		    double minosdirectionStart1[3]={  trk_dcosx_minos,  trk_dcosy_minos,  trk_dcosz_minos };
		    
		    double minosStart1[3]={ trk_vtxx_minos,  trk_vtxy_minos,  trk_vtxz_minos };
		    double xdiff1,ydiff1,rdiff1,totaldiff1,thetadiff1;
		    
		    
		    project(lardirectionEnd1, larEnd1, minosdirectionStart1,minosStart1,xdiff1,ydiff1,rdiff1,totaldiff1,thetadiff1);
		    
		    if(rdiff1>RDIFF_CUT||thetadiff1>THETADIFF_CUT){
		      //std::cout << "Event that is cut out by rdiff and thetadiff cut = " << event << " with mu KE = " << lep_mom_truth << std::endl;
		      continue;
		    }
		    
		    histeff55->Fill(lep_mom_truth);
		    histeff66->Fill(lep_mom_truth);

		    histeffT55->Fill((180./3.14159)*TMath::ACos(dcosztruth));
		    histeffT66->Fill((180./3.14159)*TMath::ACos(dcosztruth));
		    if(lep_mom_truth>1)
		      g++; // number of cc nu-mu events with a reco muon starting in FV + reco vertex in FV + muon exits the TPC + minos tracks>0 + atleast 1 -vely charged MINOS track + matched with correct sign + pass rdiff|rtheta cut
		  }
		}
	      }
	    }
	  }
	}
      }

/////////////////////////////////////////////////
      
      //Now using the same cuts as used in CCInclusive Analysis

      if(
	 nuvtxx_truth>0.+FIDX_CUT 
	 && nuvtxx_truth<47.-FIDX_CUT 
	 && nuvtxy_truth>-20.+FIDY_CUT 
	 && nuvtxy_truth<20.-FIDY_CUT 
	 && nuvtxz_truth>0.+FIDZup_CUT 
	 && nuvtxz_truth<90.-FIDZdown_CUT 
	 && ccnc_truth==0 && nuPDG_truth==14
	 && test_charge_minos==-1
	 ) {
	
	histeff1111D->Fill(lep_mom_truth);
	histeffT1111D->Fill((180./3.14159)*TMath::ACos(dcosztruth));
	if(lep_mom_truth>1)
	  j++;
      }

      if(
	 vtxx_reco>0.+FIDX_CUT &&  trackstart_x_reco>0.+FIDX_CUT 
	 &&  vtxx_reco<47.-FIDX_CUT &&  trackstart_x_reco<47.-FIDX_CUT  
	 &&  vtxy_reco>-20.+FIDY_CUT &&  trackstart_y_reco>-20.+FIDY_CUT 
	 &&  vtxy_reco<20.-FIDY_CUT  &&  trackstart_y_reco<20.-FIDY_CUT
	 &&  vtxz_reco>0.+FIDZup_CUT &&  trackstart_z_reco>0.+FIDZup_CUT 
	 &&  vtxz_reco<90.-FIDZdown_CUT &&  trackstart_z_reco<90.-FIDZdown_CUT ){
	
	//obvious cuts
	if( nmatched_reco==1&& trk_charge_minos<0&& trk_mom_minos>0){

	  double lardirectionEnd[3]={ trackexit_dcosx_reco,  trackexit_dcosy_reco,  trackexit_dcosz_reco};
	  double larEnd[3]= {  trackexit_x_reco,  trackexit_y_reco,  trackexit_z_reco};
	  double minosdirectionStart[3]={  trk_dcosx_minos,  trk_dcosy_minos,  trk_dcosz_minos };
	    
	  double minosStart[3]={ trk_vtxx_minos,  trk_vtxy_minos,  trk_vtxz_minos };
	  double xdiff,ydiff,rdiff,totaldiff,thetadiff;
	    
	  
	  project(lardirectionEnd, larEnd, minosdirectionStart,minosStart,xdiff,ydiff,rdiff,totaldiff,thetadiff);
	  
	  if(rdiff>RDIFF_CUT||thetadiff>THETADIFF_CUT)
	    continue;
	  
	  if(ccnc_truth==0 && nuPDG_truth==14){
	      histeff2222->Fill(lep_mom_truth);
	      histeffT2222->Fill((180./3.14159)*TMath::ACos(dcosztruth));
	      if(lep_mom_truth>1)
	      k++;
	  }
	}	
      }      
   }
   
   //std::cout << "a = " << a << ", b = " << b <<", c = " << c << ", d = " << d << ", e = " << e <<", f = " << f << ", g = " << g << ", j = " << j << ", k = " << k << ", true = " << true << std::endl;

   new TCanvas;
   histeff00->Divide(histeff00, histeff0D, 1., 1., "B");
   histeff00->SetTitle("CC #nu_{#mu} Muon Track Reconstruction Probability in ArgoNeuT (FV cut on Reco Muon Track Start) ");
   histeff00->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeff00->SetYTitle("probability");
   histeff00->SetMaximum(1.01);
   histeff00->SetMinimum(0);
   histeff00->Draw();
      
   new TCanvas;
   histeff11->Divide(histeff11, histeff1D, 1., 1., "B");
   histeff11->SetTitle("CC #nu_{#mu} Vertex Reconstruction Probability (after Muon Track Reconstruction)");
   histeff11->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeff11->SetYTitle("probability");
   histeff11->SetMaximum(1.01);
   histeff11->SetMinimum(0);
   histeff11->Draw();
      
   new TCanvas;
   histeff22->Divide(histeff22, histeff2D, 1., 1., "B");
   histeff22->SetTitle("CC #nu_{#mu} Muon Exiting ArgoNeuT Probability (after Muon Track Reconstruction + Vertex Reconstruction)");
   histeff22->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeff22->SetYTitle("probability");
   histeff22->SetMaximum(1.01);
   histeff22->SetMinimum(0.);
   histeff22->Draw();
    
   new TCanvas;
   histeff33->Divide(histeff33, histeff3D, 1., 1., "B");
   histeff33->SetTitle("CC #nu_{#mu} MINOS reconstruction Probability (after Muon Track Reconstruction + Vertex Reconstruction + Muon Exiting ArgoNeuT)");
   histeff33->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeff33->SetYTitle("probability");
   histeff33->SetMaximum(1.01);
   histeff33->SetMinimum(0.);
   histeff33->Draw();
    
   new TCanvas;
   histeffXX->Divide(histeffXX, histeff4D, 1., 1., "B");
   histeffXX->SetTitle("CC #nu_{#mu} MINOS Reconstruction Probability WITH correct charge (after Muon Track Reconstruction + Vertex Reconstruction + Muon Exiting ArgoNeuT + >0 MINOS Tracks)");
   histeffXX->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeffXX->SetYTitle("probability");
   histeffXX->SetMaximum(1.01);
   histeffXX->SetMinimum(0.);
   histeffXX->Draw();
    
   new TCanvas;
   histeff44->Divide(histeff44, histeffXD, 1., 1., "B");
   histeff44->SetTitle("CC #nu_{#mu} Matching Probability (after Muon Track Reconstruction + Vertex Reconstruction + Muon Exiting ArgoNeuT + >0 MINOS Tracks + correct sign)");
   histeff44->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeff44->SetYTitle("probability");
   histeff44->SetMaximum(1.01);
   histeff44->SetMinimum(0.);
   histeff44->Draw();
    
 //new TCanvas;
   histeff55->Divide(histeff55, histeffXD, 1., 1., "B");
   histeff55->SetTitle("CC #nu_{#mu} Matching Probability (after Muon Track Reconstruction + Vertex Reconstruction + Muon Exiting ArgoNeuT + >0 MINOS Tracks + correct sign)");
   histeff55->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeff55->SetYTitle("probability");
   histeff55->SetMaximum(1.01);
   histeff55->SetMinimum(0.);
   histeff55->SetLineColor(2); 
   histeff55->Draw("SAMES");
    
   new TCanvas;
   histeff66->Divide(histeff66, histeff0D, 1., 1., "B");
   histeff66->SetTitle("CC #nu_{#mu} Muon Total reconstruction + matching Probability (after all the previous cuts)");
   histeff66->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeff66->SetYTitle("probability");
   histeff66->SetMaximum(1.01);
   histeff66->SetMinimum(0.);
   histeff66->Draw();
 
   new TCanvas;
   histeff2222->Divide(histeff2222, histeff1111D, 1, 1, "B");
   histeff2222->SetTitle("CC #nu_{#mu} Muon ArgoNeuT+matching reconstruction probability (after cuts)");
   histeff2222->SetXTitle("P_{#mu, true} (GeV/c) ");
   histeff2222->SetYTitle("probability");
   histeff2222->SetMaximum(1.01);
   histeff2222->SetMinimum(0);
   histeff2222->Draw();
   

 /////for muon angle
   new TCanvas;
   histeffT00->Divide(histeffT00, histeffT0D, 1., 1., "B");
   histeffT00->SetTitle("CC #nu_{#mu} Muon Track Reconstruction Probability in ArgoNeuT (FV cut on Reco Muon Track Start) ");
   histeffT00->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffT00->SetYTitle("probability");
   histeffT00->SetMaximum(1.01);
   histeffT00->SetMinimum(0);
   histeffT00->Draw();
      
   new TCanvas;
   histeffT11->Divide(histeffT11, histeffT1D, 1., 1., "B");
   histeffT11->SetTitle("CC #nu_{#mu} Vertex Reconstruction Probability (after Muon Track Reconstruction)");
   histeffT11->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffT11->SetYTitle("probability");
   histeffT11->SetMaximum(1.01);
   histeffT11->SetMinimum(0);
   histeffT11->Draw();
      
   new TCanvas;
   histeffT22->Divide(histeffT22, histeffT2D, 1., 1., "B");
   histeffT22->SetTitle("CC #nu_{#mu} Muon Exiting ArgoNeuT Probability (after Muon Track Reconstruction + Vertex Reconstruction)");
   histeffT22->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffT22->SetYTitle("probability");
   histeffT22->SetMaximum(1.01);
   histeffT22->SetMinimum(0.);
   histeffT22->Draw();
    
   new TCanvas;
   histeffT33->Divide(histeffT33, histeffT3D, 1., 1., "B");
   histeffT33->SetTitle("CC #nu_{#mu} MINOS reconstruction Probability (after Muon Track Reconstruction + Vertex Reconstruction + Muon Exiting ArgoNeuT)");
   histeffT33->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffT33->SetYTitle("probability");
   histeffT33->SetMaximum(1.01);
   histeffT33->SetMinimum(0.);
   histeffT33->Draw();
    
   new TCanvas;
   histeffTXX->Divide(histeffTXX, histeffT4D, 1., 1., "B");
   histeffTXX->SetTitle("CC #nu_{#mu} MINOS Reconstruction Probability WITH correct charge (after Muon Track Reconstruction + Vertex Reconstruction + Muon Exiting ArgoNeuT + >0 MINOS Tracks)");
   histeffTXX->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffTXX->SetYTitle("probability");
   histeffTXX->SetMaximum(1.01);
   histeffTXX->SetMinimum(0.);
   histeffTXX->Draw();
    
   new TCanvas;
   histeffT44->Divide(histeffT44, histeffTXD, 1., 1., "B");
   histeffT44->SetTitle("CC #nu_{#mu} Matching Probability (after Muon Track Reconstruction + Vertex Reconstruction + Muon Exiting ArgoNeuT + >0 MINOS Tracks + correct sign)");
   histeffT44->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffT44->SetYTitle("probability");
   histeffT44->SetMaximum(1.01);
   histeffT44->SetMinimum(0.);
   histeffT44->Draw();
    
 //new TCanvas;
   histeffT55->Divide(histeffT55, histeffTXD, 1., 1., "B");
   histeffT55->SetTitle("CC #nu_{#mu} Matching Probability (after Muon Track Reconstruction + Vertex Reconstruction + Muon Exiting ArgoNeuT + >0 MINOS Tracks + correct sign)");
   histeffT55->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffT55->SetYTitle("probability");
   histeffT55->SetMaximum(1.01);
   histeffT55->SetMinimum(0.);
   histeffT55->SetLineColor(2); 
   histeffT55->Draw("SAMES");
    
   new TCanvas;
   histeffT66->Divide(histeffT66, histeffT0D, 1., 1., "B");
   histeffT66->SetTitle("CC #nu_{#mu} Muon Total reconstruction + matching Probability (after all the previous cuts)");
   histeffT66->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffT66->SetYTitle("probability");
   histeffT66->SetMaximum(1.01);
   histeffT66->SetMinimum(0.);
   histeffT66->Draw();
 
   new TCanvas;
   histeffT2222->Divide(histeffT2222, histeffT1111D, 1, 1, "B");
   histeffT2222->SetTitle("CC #nu_{#mu} Muon ArgoNeuT+matching reconstruction probability (after cuts)");
   histeffT2222->SetXTitle("#theta_{#mu, true} (degrees) ");
   histeffT2222->SetYTitle("probability");
   histeffT2222->SetMaximum(1.01);
   histeffT2222->SetMinimum(0);
   histeffT2222->Draw();

}


void project(double lardirectionEnd[3], double larEnd[3], double minosdirectionStart[3], double minosStart[3], double &xpred, double &ypred, double &rdiff, double &totaldiff, double &thetadiff)
{
  double D=(90*0.5)+(42.4*2.54)-5.588; //distance from the front (upstream) of the TPC to the 1st Minos plane 
  double x_offset=117.4; // previously 116.9;
  double y_offset=19.3; // previously  20.28;
  
  double dz = D - larEnd[2]+(100.0 * minosStart[2]);//z-difference between end of T962 track and
  //begin of MINOS track...in centimeters
  double l = dz/(lardirectionEnd[2]);//3-d distance between end of T962 track and begin of MINOS track
  double x_pred = l*lardirectionEnd[0]+larEnd[0];//predicted x-pos. of T962 track at z-position equal to                                                       //start of MINOS track
  double y_pred = l*lardirectionEnd[1]+larEnd[1];//predicted y-pos. of T962 track at z-position equal to                                                       //start of MINOS track
  
  double dx = 100.0*minosStart[0]- x_offset - x_pred;
  double dy = 100.0*minosStart[1] + y_offset - y_pred;
  
  
  xpred = dx;
  ypred = dy;
  rdiff = sqrt(dx*dx + dy*dy);
  //totaldiff is a measure of the agreement between the ArgoNeuT projected track and the MINOS track based on radial distance and angle. totaldiff= rdiff/cos(theta)  
  totaldiff=fabs(rdiff/((lardirectionEnd[0]*minosdirectionStart[0])+(lardirectionEnd[1]*minosdirectionStart[1])+(lardirectionEnd[2]*minosdirectionStart[2])));
  thetadiff =TMath::ACos((lardirectionEnd[0]*minosdirectionStart[0])+(lardirectionEnd[1]*minosdirectionStart[1])+(lardirectionEnd[2]*minosdirectionStart[2])); 
  
  return;
}
   
