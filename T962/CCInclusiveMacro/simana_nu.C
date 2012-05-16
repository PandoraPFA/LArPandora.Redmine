#define simkinana_cxx
#include "simkinana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "cutvariables.h"
#include "plot.C"

// main analysis loop
void simkinana::Loop()
{
  // Histogram Style Settings 
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1,0); 
  gStyle->SetLineWidth(2); 
  gStyle->SetHistLineWidth(3);
  gStyle->SetOptFit(kTRUE);
  gStyle->SetOptStat(1000001110);
  gStyle->SetOptFit(0011);
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(.032,"Z");
  gStyle->SetStatFontSize(.02);
  gStyle->SetTitleFontSize(.044);
  TH1::AddDirectory(false);
  set_plot_style();

	//////////////////////////////////////////////////////////////////
  // initialize histograms
	//////////////////////////////////////////////////////////////////

	TH1D*  evt_truF_rec0_pNU_energy   = new TH1D("", "",60,0,30);
	TH1D*  evt_truF_rec0_pANU_energy   = new TH1D("", "",60,0,30);
	
	// true fiducial volume
  double evt_truF_rec0_p0_count	   = 0;
  TH1D*  evt_truF_rec0_p0_theta    = new TH1D("", "",18,0.,36);     evt_truF_rec0_p0_theta->Sumw2();
  TH1D*  evt_truF_rec0_p0_momentum = new TH1D("", "",20,0.,25);  evt_truF_rec0_p0_momentum->Sumw2();
	TH1D*  evt_truF_rec0_p0_energy   = new TH1D("", "",60,0,30);
	// true fiducial volume, true CC
  double evt_truF_rec0_pCC_count       = 0;
  double evt_truF_rec0_pCC_count_mom   = 0;
  double evt_truF_rec0_pCC_count_theta = 0;
  TH1D*  evt_truF_rec0_pCC_theta       = new TH1D("", "",18,0.,36); evt_truF_rec0_pCC_theta->Sumw2();
  TH1D*  evt_truF_rec0_pCC_momentum    = new TH1D("", "",20,0.,25); evt_truF_rec0_pCC_momentum->Sumw2();
	TH1D*  evt_truF_rec0_pCC_energy      = new TH1D("", "",60,0,30);
  // ccqe
  double evt_truF_rec0_pCCQE_count    = 0;
  TH1D*  evt_truF_rec0_pCCQE_theta    = new TH1D("", "",18,0.,36); evt_truF_rec0_pCCQE_theta->Sumw2();
  TH1D*  evt_truF_rec0_pCCQE_momentum = new TH1D("", "",20,0.,25); evt_truF_rec0_pCCQE_momentum->Sumw2();
  TH1D*  evt_truF_rec0_pCCQE_energy   = new TH1D("", "",60,0,30);
	// ccres
  double evt_truF_rec0_pCCRES_count    = 0;
  TH1D*  evt_truF_rec0_pCCRES_theta    = new TH1D("", "",18,0.,36); evt_truF_rec0_pCCRES_theta->Sumw2();
  TH1D*  evt_truF_rec0_pCCRES_momentum = new TH1D("", "",20,0.,25); evt_truF_rec0_pCCRES_momentum->Sumw2();
  TH1D*  evt_truF_rec0_pCCRES_energy   = new TH1D("", "",60,0,30);
	// ccdis
  double evt_truF_rec0_pCCDIS_count    = 0;
  TH1D*  evt_truF_rec0_pCCDIS_theta    = new TH1D("", "",18,0.,36); evt_truF_rec0_pCCDIS_theta->Sumw2();
  TH1D*  evt_truF_rec0_pCCDIS_momentum = new TH1D("", "",20,0.,25); evt_truF_rec0_pCCDIS_momentum->Sumw2();
  TH1D*  evt_truF_rec0_pCCDIS_energy   = new TH1D("", "",60,0,30);
	// ccdis
  double evt_truF_rec0_pCCOTHER_count    = 0;
  TH1D*  evt_truF_rec0_pCCOTHER_theta    = new TH1D("", "",18,0.,36); evt_truF_rec0_pCCOTHER_theta->Sumw2();
  TH1D*  evt_truF_rec0_pCCOTHER_momentum = new TH1D("", "",20,0.,25); evt_truF_rec0_pCCOTHER_momentum->Sumw2();
  TH1D*  evt_truF_rec0_pCCOTHER_energy   = new TH1D("", "",60,0,30);
	// BKG
  double evt_truF_rec0_pBKG_count    = 0;
  TH1D*  evt_truF_rec0_pBKG_theta    = new TH1D("", "",18,0.,36); evt_truF_rec0_pBKG_theta->Sumw2();
  TH1D*  evt_truF_rec0_pBKG_momentum = new TH1D("", "",20,0.,25); evt_truF_rec0_pBKG_momentum->Sumw2();
  TH1D*  evt_truF_rec0_pBKG_energy   = new TH1D("", "",60,0,30);
	// nc
  double evt_truF_rec0_pNC_count    = 0;
  TH1D*  evt_truF_rec0_pNC_theta    = new TH1D("", "",18,0.,36); evt_truF_rec0_pNC_theta->Sumw2();
  TH1D*  evt_truF_rec0_pNC_momentum = new TH1D("", "",20,0.,25); evt_truF_rec0_pNC_momentum->Sumw2();
  TH1D*  evt_truF_rec0_pNC_energy   = new TH1D("", "",60,0,30);
	// WS
  double evt_truF_rec0_pWS_count    = 0;
  TH1D*  evt_truF_rec0_pWS_theta    = new TH1D("", "",18,0.,36); evt_truF_rec0_pWS_theta->Sumw2();
  TH1D*  evt_truF_rec0_pWS_momentum = new TH1D("", "",20,0.,25); evt_truF_rec0_pWS_momentum->Sumw2();
  TH1D*  evt_truF_rec0_pWS_energy   = new TH1D("", "",60,0,30);
	// other
  double evt_truF_rec0_pBOTHER_count    = 0;
  TH1D*  evt_truF_rec0_pBOTHER_theta    = new TH1D("", "",18,0.,36); evt_truF_rec0_pBOTHER_theta->Sumw2();
  TH1D*  evt_truF_rec0_pBOTHER_momentum = new TH1D("", "",20,0.,25); evt_truF_rec0_pBOTHER_momentum->Sumw2();
  TH1D*  evt_truF_rec0_pBOTHER_energy   = new TH1D("", "",60,0,30);

  // true fiducial volume, true neutrino, MINOS match, MINOS charge
  double evt_truFQ_rec0_pCC_count			= 0;
  TH1D* evt_truFQ_rec0_p0_theta = new TH1D("", "",18,0.,36);
  evt_truFQ_rec0_p0_theta->Sumw2();
  TH1D* evt_truFQ_rec0_p0_momentum = new TH1D("", "",20,0.,25);
  evt_truFQ_rec0_p0_momentum->Sumw2();
  TH1D* evt_truFQ_rec0_p0_energy = new TH1D("", "",60,0,30);

  double evt_truFQ_rec0_pCCQE_count		= 0;
  double evt_truFQ_rec0_pCCDIS_count		= 0;
  double evt_truFQ_rec0_pCCRES_count		= 0;
  double evt_truFQ_rec0_pCCOTHER_count		= 0;


  //////////////////////////////////////////////////////////////////
	// reconstruction
  //////////////////////////////////////////////////////////////////

	// reco fiducial volume
  double evt_tru0_recF_p0_count    = 0;
  TH1D*  evt_tru0_recF_p0_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_p0_theta->Sumw2();
  TH1D*  evt_tru0_recF_p0_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_p0_momentum->Sumw2();
  TH1D*  evt_tru0_recF_p0_energy   = new TH1D("", "",60,0,30);

	// reco fiducial volume, CC
  double evt_tru0_recF_pCC_count    = 0;
  TH1D*  evt_tru0_recF_pCC_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pCC_theta->Sumw2();
  TH1D*  evt_tru0_recF_pCC_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pCC_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pCC_energy   = new TH1D("", "",60,0,30);
  // ccqe
  double evt_tru0_recF_pCCQE_count    = 0;
  TH1D*  evt_tru0_recF_pCCQE_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pCCQE_theta->Sumw2();
  TH1D*  evt_tru0_recF_pCCQE_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pCCQE_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pCCQE_energy   = new TH1D("", "",60,0,30);
	// ccres
  double evt_tru0_recF_pCCRES_count    = 0;
  TH1D*  evt_tru0_recF_pCCRES_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pCCRES_theta->Sumw2();
  TH1D*  evt_tru0_recF_pCCRES_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pCCRES_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pCCRES_energy   = new TH1D("", "",60,0,30);
	// ccdis
  double evt_tru0_recF_pCCDIS_count    = 0;
  TH1D*  evt_tru0_recF_pCCDIS_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pCCDIS_theta->Sumw2();
  TH1D*  evt_tru0_recF_pCCDIS_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pCCDIS_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pCCDIS_energy   = new TH1D("", "",60,0,30);
	// ccdis
  double evt_tru0_recF_pCCOTHER_count    = 0;
  TH1D*  evt_tru0_recF_pCCOTHER_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pCCOTHER_theta->Sumw2();
  TH1D*  evt_tru0_recF_pCCOTHER_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pCCOTHER_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pCCOTHER_energy   = new TH1D("", "",60,0,30);
	// BKG
  double evt_tru0_recF_pBKG_count    = 0;
  TH1D*  evt_tru0_recF_pBKG_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pBKG_theta->Sumw2();
  TH1D*  evt_tru0_recF_pBKG_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pBKG_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pBKG_energy   = new TH1D("", "",60,0,30);
	// nc
  double evt_tru0_recF_pNC_count    = 0;
  TH1D*  evt_tru0_recF_pNC_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pNC_theta->Sumw2();
  TH1D*  evt_tru0_recF_pNC_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pNC_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pNC_energy   = new TH1D("", "",60,0,30);
	// WS
  double evt_tru0_recF_pWS_count    = 0;
  TH1D*  evt_tru0_recF_pWS_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pWS_theta->Sumw2();
  TH1D*  evt_tru0_recF_pWS_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pWS_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pWS_energy   = new TH1D("", "",60,0,30);
	// other
  double evt_tru0_recF_pBOTHER_count    = 0;
  TH1D*  evt_tru0_recF_pBOTHER_theta    = new TH1D("", "",18,0.,36); evt_tru0_recF_pBOTHER_theta->Sumw2();
  TH1D*  evt_tru0_recF_pBOTHER_momentum = new TH1D("", "",20,0.,25); evt_tru0_recF_pBOTHER_momentum->Sumw2();
  TH1D*  evt_tru0_recF_pBOTHER_energy   = new TH1D("", "",60,0,30);

	// reco fiducial volume, true fiducial volume
  double evt_truF_recF_p0_count    = 0;
  TH1D*  evt_truF_recF_p0_theta    = new TH1D("", "",18,0.,36); evt_truF_recF_p0_theta->Sumw2();
  TH1D*  evt_truF_recF_p0_momentum = new TH1D("", "",20,0.,25); evt_truF_recF_p0_momentum->Sumw2();
  TH1D*  evt_truF_recF_p0_energy   = new TH1D("", "",60,0,30);

	// reco fiducial volume, MINOS match
  double evt_tru0_recFM_p0_count    = 0;
  TH1D*  evt_tru0_recFM_p0_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_p0_theta->Sumw2();
  TH1D*  evt_tru0_recFM_p0_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_p0_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_p0_energy   = new TH1D("", "",60,0,30);
	// reco: fiducial volume, MINOS match, CC
  double evt_tru0_recFM_pCC_count    = 0;
  TH1D*  evt_tru0_recFM_pCC_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pCC_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pCC_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pCC_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pCC_energy   = new TH1D("", "",60,0,30);
  // ccqe
  double evt_tru0_recFM_pCCQE_count    = 0;
  TH1D*  evt_tru0_recFM_pCCQE_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pCCQE_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pCCQE_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pCCQE_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pCCQE_energy   = new TH1D("", "",60,0,30);
	// ccres
  double evt_tru0_recFM_pCCRES_count    = 0;
  TH1D*  evt_tru0_recFM_pCCRES_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pCCRES_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pCCRES_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pCCRES_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pCCRES_energy   = new TH1D("", "",60,0,30);
	// ccdis
  double evt_tru0_recFM_pCCDIS_count    = 0;
  TH1D*  evt_tru0_recFM_pCCDIS_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pCCDIS_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pCCDIS_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pCCDIS_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pCCDIS_energy   = new TH1D("", "",60,0,30);
	// ccdis
  double evt_tru0_recFM_pCCOTHER_count    = 0;
  TH1D*  evt_tru0_recFM_pCCOTHER_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pCCOTHER_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pCCOTHER_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pCCOTHER_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pCCOTHER_energy   = new TH1D("", "",60,0,30);
	// BKG
  double evt_tru0_recFM_pBKG_count    = 0;
  TH1D*  evt_tru0_recFM_pBKG_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pBKG_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pBKG_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pBKG_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pBKG_energy   = new TH1D("", "",60,0,30);
	// nc
  double evt_tru0_recFM_pNC_count    = 0;
  TH1D*  evt_tru0_recFM_pNC_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pNC_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pNC_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pNC_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pNC_energy   = new TH1D("", "",60,0,30);
	// WS
  double evt_tru0_recFM_pWS_count    = 0;
  TH1D*  evt_tru0_recFM_pWS_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pWS_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pWS_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pWS_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pWS_energy   = new TH1D("", "",60,0,30);
	// other
  double evt_tru0_recFM_pBOTHER_count    = 0;
  TH1D*  evt_tru0_recFM_pBOTHER_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFM_pBOTHER_theta->Sumw2();
  TH1D*  evt_tru0_recFM_pBOTHER_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFM_pBOTHER_momentum->Sumw2();
  TH1D*  evt_tru0_recFM_pBOTHER_energy   = new TH1D("", "",60,0,30);

	// reco fiducial volume, MINOS match, MINOS charge 
  double evt_tru0_recFMQ_p0_count		 = 0;
  TH1D*  evt_tru0_recFMQ_p0_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_p0_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_p0_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_p0_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_p0_energy   = new TH1D("", "",60,0,30);
 
  // cc-inclusive
  double evt_tru0_recFMQ_pCC_count	     = 0;  
  double evt_tru0_recFMQ_pCC_count_mom   = 0;
  double evt_tru0_recFMQ_pCC_count_theta = 0;  
  TH1D*  evt_tru0_recFMQ_pCC_theta       = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pCC_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCC_momentum    = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pCC_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCC_energy      = new TH1D("", "",60,0,30);
  // ccqe
  double evt_tru0_recFMQ_pCCQE_count    = 0;
  TH1D*  evt_tru0_recFMQ_pCCQE_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pCCQE_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCCQE_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pCCQE_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCCQE_energy   = new TH1D("", "",60,0,30);
	// ccres
  double evt_tru0_recFMQ_pCCRES_count    = 0;
  TH1D*  evt_tru0_recFMQ_pCCRES_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pCCRES_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCCRES_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pCCRES_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCCRES_energy   = new TH1D("", "",60,0,30);
	// ccdis
  double evt_tru0_recFMQ_pCCDIS_count    = 0;
  TH1D*  evt_tru0_recFMQ_pCCDIS_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pCCDIS_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCCDIS_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pCCDIS_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCCDIS_energy   = new TH1D("", "",60,0,30);
	// ccdis
  double evt_tru0_recFMQ_pCCOTHER_count    = 0;
  TH1D*  evt_tru0_recFMQ_pCCOTHER_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pCCOTHER_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCCOTHER_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pCCOTHER_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pCCOTHER_energy   = new TH1D("", "",60,0,30);
	// BKG
  double evt_tru0_recFMQ_pBKG_count    = 0;
  TH1D*  evt_tru0_recFMQ_pBKG_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pBKG_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pBKG_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pBKG_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pBKG_energy   = new TH1D("", "",60,0,30);
	// nc
  double evt_tru0_recFMQ_pNC_count    = 0;
  TH1D*  evt_tru0_recFMQ_pNC_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pNC_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pNC_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pNC_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pNC_energy   = new TH1D("", "",60,0,30);
	// WS
  double evt_tru0_recFMQ_pWS_count    = 0;
  TH1D*  evt_tru0_recFMQ_pWS_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pWS_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pWS_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pWS_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pWS_energy   = new TH1D("", "",60,0,30);
	// other
  double evt_tru0_recFMQ_pBOTHER_count    = 0;
  TH1D*  evt_tru0_recFMQ_pBOTHER_theta    = new TH1D("", "",18,0.,36); evt_tru0_recFMQ_pBOTHER_theta->Sumw2();
  TH1D*  evt_tru0_recFMQ_pBOTHER_momentum = new TH1D("", "",20,0.,25); evt_tru0_recFMQ_pBOTHER_momentum->Sumw2();
  TH1D*  evt_tru0_recFMQ_pBOTHER_energy   = new TH1D("", "",60,0,30);

	// Open AnalysisTree (use chain for larger files---for breaking into parts)
	TChain ch("KinTree");
//	ch.Add("20120119_simkinreco_kinga.root/ccqeanalysistree/anatree"); 
	ch.Add("20120307_mcrecokalman_nu.root/analysistree/anatree"); 
  simkinana aEvent(KinTree);
  Long64_t nentries = KinTree->GetEntries();
  TFile *outfile = new TFile("out.root","RECREATE");
  cout << "Number of Entries " << nentries << endl;

  // initialize event counts
  float theta												= 0.;
  double totpot;

	// read in systematics profile
  TFile *fi=new TFile("profsyst.root","READ");
  TH1F *profthetasyst2 = (TH1F*)fi->Get("histthetasystfinal");
  TH1F *profmomsyst2 = (TH1F*)fi->Get("histmomsystfinal");   
  double thetasyst[30]	= { profthetasyst2->GetBinContent(1),  profthetasyst2->GetBinContent(2),  profthetasyst2->GetBinContent(3), 
  													profthetasyst2->GetBinContent(4),  profthetasyst2->GetBinContent(5),  profthetasyst2->GetBinContent(6), 
  													profthetasyst2->GetBinContent(7),  profthetasyst2->GetBinContent(8),  profthetasyst2->GetBinContent(9), 
  													profthetasyst2->GetBinContent(10), profthetasyst2->GetBinContent(11), profthetasyst2->GetBinContent(12), 
  													profthetasyst2->GetBinContent(13), profthetasyst2->GetBinContent(14), profthetasyst2->GetBinContent(15),
  													profthetasyst2->GetBinContent(16), profthetasyst2->GetBinContent(17), profthetasyst2->GetBinContent(18),
  													profthetasyst2->GetBinContent(19), profthetasyst2->GetBinContent(20), profthetasyst2->GetBinContent(21),
  													profthetasyst2->GetBinContent(22), profthetasyst2->GetBinContent(23), profthetasyst2->GetBinContent(24),
  													profthetasyst2->GetBinContent(25), profthetasyst2->GetBinContent(26), profthetasyst2->GetBinContent(27),
  													profthetasyst2->GetBinContent(28), profthetasyst2->GetBinContent(29), profthetasyst2->GetBinContent(30)};

  double momsyst[30]		= { profmomsyst2->GetBinContent(1),  profmomsyst2->GetBinContent(2),  profmomsyst2->GetBinContent(3), 
  													profmomsyst2->GetBinContent(4),  profmomsyst2->GetBinContent(5),  profmomsyst2->GetBinContent(6), 
  													profmomsyst2->GetBinContent(7),  profmomsyst2->GetBinContent(8),  profmomsyst2->GetBinContent(9), 
  													profmomsyst2->GetBinContent(10), profmomsyst2->GetBinContent(11), profmomsyst2->GetBinContent(12), 
  													profmomsyst2->GetBinContent(13), profmomsyst2->GetBinContent(14), profmomsyst2->GetBinContent(15),
  													profmomsyst2->GetBinContent(16), profmomsyst2->GetBinContent(17), profmomsyst2->GetBinContent(18),
  													profmomsyst2->GetBinContent(19), profmomsyst2->GetBinContent(20), profmomsyst2->GetBinContent(21),
  													profmomsyst2->GetBinContent(22), profmomsyst2->GetBinContent(23), profmomsyst2->GetBinContent(24),
  													profmomsyst2->GetBinContent(25), profmomsyst2->GetBinContent(26), profmomsyst2->GetBinContent(27),
  													profmomsyst2->GetBinContent(28), profmomsyst2->GetBinContent(29), profmomsyst2->GetBinContent(30)};

	// read in energy loss as function of distance
  TFile *f = new TFile("energylost.root","READ");
  TH1F *energylost = (TH1F*)f->Get("energylost");
  double energyloss[14]	= { energylost->GetBinContent(1),  energylost->GetBinContent(2),  energylost->GetBinContent(3), 
  													energylost->GetBinContent(4),  energylost->GetBinContent(5),  energylost->GetBinContent(6), 
  													energylost->GetBinContent(7),  energylost->GetBinContent(8),  energylost->GetBinContent(9), 
  													energylost->GetBinContent(10), energylost->GetBinContent(11), energylost->GetBinContent(12), 
  													energylost->GetBinContent(13), energylost->GetBinContent(14)};

	// read in numu_numode_final.root
	// do I want to do this for anti-neutrino background simulation as well??
  TFile *f = new TFile("numu_numode_final.root","READ");
  TH1F *reweight = (TH1F*)f->Get("histdiv");
	double rweight=0.;  
  double weight[14]	= { reweight->GetBinContent(1),  reweight->GetBinContent(2),  reweight->GetBinContent(3), 
												reweight->GetBinContent(4),  reweight->GetBinContent(5),  reweight->GetBinContent(6), 
												reweight->GetBinContent(7),  reweight->GetBinContent(8),  reweight->GetBinContent(9), 
												reweight->GetBinContent(10), reweight->GetBinContent(11), reweight->GetBinContent(12), 
												reweight->GetBinContent(13), reweight->GetBinContent(14)};
  
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// LOOP OVER ALL EVENTS
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int event = 0 ; event < nentries ; event++) {
	
		if ( event % 10000 == 0 ) cout << event << " / " << nentries << endl;
		
		aEvent.GetEntry(event);  
		totpot += aEvent.pot / 100.;		
		
		// reconstructed-Y direction cosign 
		double dcosy			= ( aEvent.trackstart_dcosy_reco * TMath::Cos(.0583497) ) 
		                      + ( aEvent.trackstart_dcosz_reco * TMath::Sin(.0583497) );

		// reconstructed-Z direction cosign 
		double dcosz			= ( aEvent.trackstart_dcosz_reco * TMath::Cos(.0583497) ) 
		                      - ( aEvent.trackstart_dcosy_reco * TMath::Sin(.0583497) );

		// true-Y direction cosign 
		double dcosytruth	= ( aEvent.lep_dcosy_truth * TMath::Cos(.0583497) ) 
		                      + ( aEvent.lep_dcosz_truth * TMath::Sin(.0583497) );

		// true-Z direction cosign		                      
		double dcosztruth	= ( aEvent.lep_dcosz_truth * TMath::Cos(.0583497) ) 
		                      - ( aEvent.lep_dcosy_truth * TMath::Sin(.0583497) );

		double rweight = Weigh( aEvent.enu_truth , weight );
		
		// iso weight p,n
		double isoweight=1.;
		if ( aEvent.hitnuc_truth == 2212 ) isoweight = 22. / 20.; // proton
		if ( aEvent.hitnuc_truth == 2112 ) isoweight = 18. / 20.; // neutron

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// true fiducial volume cut
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(			aEvent.nuvtxx_truth >   0. + FIDX_CUT 
				&& 	aEvent.nuvtxx_truth <  47. - FIDX_CUT 
				&& 	aEvent.nuvtxy_truth > -20. + FIDY_CUT 
				&& 	aEvent.nuvtxy_truth <  20. - FIDY_CUT 
				&& 	aEvent.nuvtxz_truth >   0. + FIDZup_CUT 
				&& 	aEvent.nuvtxz_truth <  90. - FIDZdown_CUT ) {
				
			evt_truF_rec0_p0_count += ( 1. * rweight );
			evt_truF_rec0_p0_energy->Fill(aEvent.enu_truth,rweight);
			evt_truF_rec0_p0_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
			evt_truF_rec0_p0_momentum->Fill(aEvent.lep_mom_truth,rweight);
			
			if ( aEvent.nuPDG_truth ==  14 ) { evt_truF_rec0_pNU_energy->Fill(aEvent.enu_truth,rweight); }
			if ( aEvent.nuPDG_truth == -14 ) { evt_truF_rec0_pANU_energy->Fill(aEvent.enu_truth,rweight); }

      if ( aEvent.ccnc_truth != 0 || aEvent.nuPDG_truth != 14 ) {
				evt_truF_rec0_pBKG_count += ( 1. * rweight );
				evt_truF_rec0_pBKG_energy->Fill(aEvent.enu_truth,rweight);
				evt_truF_rec0_pBKG_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
				evt_truF_rec0_pBKG_momentum->Fill(aEvent.lep_mom_truth,rweight);

				if ( aEvent.ccnc_truth == 1 ) {				
					evt_truF_rec0_pNC_count += ( 1. * rweight );
					evt_truF_rec0_pNC_energy->Fill(aEvent.enu_truth,rweight);
					evt_truF_rec0_pNC_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truF_rec0_pNC_momentum->Fill(aEvent.lep_mom_truth,rweight);				
				} else if ( aEvent.nuPDG_truth == -14 ) {
					evt_truF_rec0_pWS_count += ( 1. * rweight );
					evt_truF_rec0_pWS_energy->Fill(aEvent.enu_truth,rweight);
					evt_truF_rec0_pWS_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truF_rec0_pWS_momentum->Fill(aEvent.lep_mom_truth,rweight);				
				} else if ( aEvent.nuPDG_truth != 14 ) {
					evt_truF_rec0_pBOTHER_count += ( 1. * rweight );
					evt_truF_rec0_pBOTHER_energy->Fill(aEvent.enu_truth,rweight);
					evt_truF_rec0_pBOTHER_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truF_rec0_pBOTHER_momentum->Fill(aEvent.lep_mom_truth,rweight);
				}

      } else if ( aEvent.ccnc_truth == 0 && aEvent.nuPDG_truth == 14 ) {
				evt_truF_rec0_pCC_count += ( 1. * rweight );
				if( aEvent.lep_mom_truth < 25. )														evt_truF_rec0_pCC_count_mom	 += ( 1. * rweight );
				if( ( 180. / 3.14159 ) * TMath::ACos( dcosztruth ) < 36. ) 	evt_truF_rec0_pCC_count_theta += ( 1. * rweight );			
				evt_truF_rec0_pCC_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
				evt_truF_rec0_pCC_momentum->Fill(aEvent.lep_mom_truth,rweight);
				evt_truF_rec0_pCC_energy->Fill(aEvent.enu_truth,rweight);
				if ( aEvent.mode_truth == 0 ) {											// ccqe
					evt_truF_rec0_pCCQE_count += ( 1. * rweight );
					evt_truF_rec0_pCCQE_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truF_rec0_pCCQE_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_truF_rec0_pCCQE_energy->Fill(aEvent.enu_truth,rweight);
				} else if ( aEvent.mode_truth == 1 ) { 							// ccres
					evt_truF_rec0_pCCRES_count += ( 1. * rweight );
					evt_truF_rec0_pCCRES_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truF_rec0_pCCRES_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_truF_rec0_pCCRES_energy->Fill(aEvent.enu_truth,rweight);
				} else if ( aEvent.mode_truth == 2 ) {							// ccdis
					evt_truF_rec0_pCCDIS_count += ( 1. * rweight );
					evt_truF_rec0_pCCDIS_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truF_rec0_pCCDIS_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_truF_rec0_pCCDIS_energy->Fill(aEvent.enu_truth,rweight);
				} else if ( aEvent.mode_truth == 3 ) {							// ccdis
					evt_truF_rec0_pCCOTHER_count += ( 1. * rweight );
					evt_truF_rec0_pCCOTHER_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truF_rec0_pCCOTHER_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_truF_rec0_pCCOTHER_energy->Fill(aEvent.enu_truth,rweight);
				}

				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				// "true" MINOS "cut"
				//////////////////////////////////////////////////////////////////////////////////////////////////////////
				// test charge is matched track's charge correct according to MINOS (0 if no match)
				// maybe a different check here would be good to step through match and THEN charge...?
				if( aEvent.test_charge_minos == -1 ) {
					evt_truFQ_rec0_pCC_count += ( 1. * rweight ); // maybe i'm wrong here: "nu mu cc events before cuts reco MINOS" is what it seems to be saying, meaning no best fit cuts?
					evt_truFQ_rec0_p0_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truFQ_rec0_p0_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_truFQ_rec0_p0_energy->Fill(aEvent.enu_truth,rweight); 
					if 			( aEvent.mode_truth == 0 ) evt_truFQ_rec0_pCCQE_count    += ( 1. * rweight );
					else if ( aEvent.mode_truth == 1 ) evt_truFQ_rec0_pCCRES_count   += ( 1. * rweight );
					else if ( aEvent.mode_truth == 2 ) evt_truFQ_rec0_pCCDIS_count   += ( 1. * rweight );
					else if ( aEvent.mode_truth == 3 ) evt_truFQ_rec0_pCCOTHER_count += ( 1. * rweight );
				}
			}
		}


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Reconstructed information
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// fiducial volume -> MINOS match, charge, forward momentum -> best fit with match
		if (		aEvent.vtxx_reco >   0. +     FIDX_CUT   && aEvent.trackstart_x_reco >   0. +     FIDX_CUT 
				&&  aEvent.vtxx_reco <  47. -     FIDX_CUT   && aEvent.trackstart_x_reco <  47. -     FIDX_CUT  
				&&  aEvent.vtxy_reco > -20. +     FIDY_CUT   && aEvent.trackstart_y_reco > -20. +     FIDY_CUT 
				&&  aEvent.vtxy_reco <  20. -     FIDY_CUT   && aEvent.trackstart_y_reco <  20. -     FIDY_CUT
				&&  aEvent.vtxz_reco >   0. +   FIDZup_CUT   && aEvent.trackstart_z_reco >   0. +   FIDZup_CUT 
				&&  aEvent.vtxz_reco <  90. - FIDZdown_CUT   && aEvent.trackstart_z_reco <  90. - FIDZdown_CUT ) {

			evt_tru0_recF_p0_count += ( 1. * rweight );
			evt_tru0_recF_p0_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
			evt_tru0_recF_p0_momentum->Fill(aEvent.lep_mom_truth,rweight);
			evt_tru0_recF_p0_energy->Fill(aEvent.enu_truth,rweight);

			if ( aEvent.ccnc_truth != 0 || aEvent.nuPDG_truth != 14 ) {
				evt_tru0_recF_pBKG_count += ( 1. * rweight );
				evt_tru0_recF_pBKG_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
				evt_tru0_recF_pBKG_momentum->Fill(aEvent.lep_mom_truth,rweight);
				evt_tru0_recF_pBKG_energy->Fill(aEvent.enu_truth,rweight);									
				if ( aEvent.ccnc_truth == 1 ) { //NC
					evt_tru0_recF_pNC_count += ( 1. * rweight );
					evt_tru0_recF_pNC_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_tru0_recF_pNC_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_tru0_recF_pNC_energy->Fill(aEvent.enu_truth,rweight);									
				}	else if ( aEvent.nuPDG_truth == -14 ) { // WS
					evt_tru0_recF_pWS_count += ( 1. * rweight ); 
					evt_tru0_recF_pWS_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_tru0_recF_pWS_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_tru0_recF_pWS_energy->Fill(aEvent.enu_truth,rweight);
				}	else if ( aEvent.nuPDG_truth != 14 ) { // not NC, WS, or CC
					evt_tru0_recF_pBOTHER_count += ( 1. * rweight ); 
					evt_tru0_recF_pBOTHER_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_tru0_recF_pBOTHER_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_tru0_recF_pBOTHER_energy->Fill(aEvent.enu_truth,rweight);
				}
			} else if ( aEvent.ccnc_truth == 0 && aEvent.nuPDG_truth == 14 ) {
				evt_tru0_recF_pCC_count += ( 1. * rweight );
				evt_tru0_recF_pCC_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
				evt_tru0_recF_pCC_momentum->Fill(aEvent.lep_mom_truth,rweight);
				evt_tru0_recF_pCC_energy->Fill(aEvent.enu_truth,rweight);
				// interaction type
				if ( aEvent.mode_truth == 0 ) {
					evt_tru0_recF_pCCQE_count += ( 1. * rweight ); 
					evt_tru0_recF_pCCQE_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_tru0_recF_pCCQE_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_tru0_recF_pCCQE_energy->Fill(aEvent.enu_truth,rweight);
				} else if ( aEvent.mode_truth == 1 ) {
					evt_tru0_recF_pCCRES_count += ( 1. * rweight );
					evt_tru0_recF_pCCRES_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_tru0_recF_pCCRES_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_tru0_recF_pCCRES_energy->Fill(aEvent.enu_truth,rweight);
				} else if ( aEvent.mode_truth == 2 ) {
					evt_tru0_recF_pCCDIS_count += ( 1. * rweight );
					evt_tru0_recF_pCCDIS_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_tru0_recF_pCCDIS_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_tru0_recF_pCCDIS_energy->Fill(aEvent.enu_truth,rweight);
				} else if ( aEvent.mode_truth == 3 ) {
					evt_tru0_recF_pCCOTHER_count += ( 1. * rweight );
					evt_tru0_recF_pCCOTHER_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_tru0_recF_pCCOTHER_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_tru0_recF_pCCOTHER_energy->Fill(aEvent.enu_truth,rweight);
				}
			}
			
			if (		aEvent.nuvtxx_truth >   0. + FIDX_CUT 
					&& 	aEvent.nuvtxx_truth <  47. - FIDX_CUT 
					&& 	aEvent.nuvtxy_truth > -20. + FIDY_CUT 
					&& 	aEvent.nuvtxy_truth <  20. - FIDY_CUT 
					&& 	aEvent.nuvtxz_truth >   0. + FIDZup_CUT 
					&& 	aEvent.nuvtxz_truth <  90. - FIDZdown_CUT ) {
				if ( aEvent.ccnc_truth == 0 && aEvent.nuPDG_truth == 14 ) {
					evt_truF_recF_p0_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_truF_recF_p0_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_truF_recF_p0_energy->Fill(aEvent.enu_truth,rweight);
				}
			}
			
			// ArgoNeuT exit to MINOS enter angular differences
			double diff_dcosx = aEvent.trackexit_dcosx_reco - aEvent.trk_dcosx_minos;
			double diff_dcosy = aEvent.trackexit_dcosy_reco - aEvent.trk_dcosy_minos;
			double diff_dcosz = aEvent.trackexit_dcosz_reco - aEvent.trk_dcosz_minos;
			
			// MINOS matching cuts: match, charge, momentum (inside is best fit to minos cut)
			if (  aEvent.nmatched_reco   == 1 
			   && aEvent.trk_mom_minos    > 0 ) {
				
				double lardirectionEnd[3]			=	{	aEvent.trackexit_dcosx_reco, 
																					aEvent.trackexit_dcosy_reco, 
																					aEvent.trackexit_dcosz_reco };

				double larEnd[3]							= { aEvent.trackexit_x_reco, 
																					aEvent.trackexit_y_reco, 
																					aEvent.trackexit_z_reco };
																			
				double minosdirectionStart[3] =	{ aEvent.trk_dcosx_minos, 
																					aEvent.trk_dcosy_minos, 
																					aEvent.trk_dcosz_minos };
				
				double minosStart[3]					= { aEvent.trk_vtxx_minos, 
																					aEvent.trk_vtxy_minos, 
																					aEvent.trk_vtxz_minos };
																					
				double xdiff, ydiff, rdiff, totaldiff, thetadiff;		
				
				project( lardirectionEnd, larEnd, minosdirectionStart, minosStart, xdiff, ydiff, rdiff, totaldiff, thetadiff );
				
				double minosvtx[3] = { ( aEvent.trk_vtxx_minos * 100. ) - 117.4, 
															 ( aEvent.trk_vtxy_minos * 100. ) +  19.3, 
															 ( aEvent.trk_vtxz_minos * 100. ) + 147.1  };
				
				double distance = sqrt(   pow( aEvent.vtxx_reco - minosvtx[0] , 2 )
				                        + pow( aEvent.vtxy_reco - minosvtx[1] , 2 )
				                        + pow( aEvent.vtxz_reco - minosvtx[2] , 2 ) );
				
				// best fit cut
				if ( rdiff < RDIFF_CUT || thetadiff < THETADIFF_CUT ) {
				
          evt_tru0_recFM_p0_count += ( 1. * rweight );
					evt_tru0_recFM_p0_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
					evt_tru0_recFM_p0_momentum->Fill(aEvent.lep_mom_truth,rweight);
					evt_tru0_recFM_p0_energy->Fill(aEvent.enu_truth,rweight);
					if ( aEvent.ccnc_truth != 0 || aEvent.nuPDG_truth != 14 ) {
						evt_tru0_recFM_pBKG_count += ( 1. * rweight );
						evt_tru0_recFM_pBKG_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
						evt_tru0_recFM_pBKG_momentum->Fill(aEvent.lep_mom_truth,rweight);
						evt_tru0_recFM_pBKG_energy->Fill(aEvent.enu_truth,rweight);									
						if ( aEvent.ccnc_truth == 1 ) { //NC
							evt_tru0_recFM_pNC_count += ( 1. * rweight );
							evt_tru0_recFM_pNC_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFM_pNC_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFM_pNC_energy->Fill(aEvent.enu_truth,rweight);									
						}	else if ( aEvent.nuPDG_truth == -14 ) { // WS
							evt_tru0_recFM_pWS_count += ( 1. * rweight ); 
							evt_tru0_recFM_pWS_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFM_pWS_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFM_pWS_energy->Fill(aEvent.enu_truth,rweight);
						}	else if ( aEvent.nuPDG_truth != 14 ) { // not NC, WS, or CC
							evt_tru0_recFM_pBOTHER_count += ( 1. * rweight ); 
							evt_tru0_recFM_pBOTHER_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFM_pBOTHER_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFM_pBOTHER_energy->Fill(aEvent.enu_truth,rweight);
						}
					} else if ( aEvent.ccnc_truth == 0 && aEvent.nuPDG_truth == 14 ) {
						evt_tru0_recFM_pCC_count += ( 1. * rweight );
						evt_tru0_recFM_pCC_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
						evt_tru0_recFM_pCC_momentum->Fill(aEvent.lep_mom_truth,rweight);
						evt_tru0_recFM_pCC_energy->Fill(aEvent.enu_truth,rweight);
						// interaction type
						if ( aEvent.mode_truth == 0 ) {
							evt_tru0_recFM_pCCQE_count += ( 1. * rweight ); 
							evt_tru0_recFM_pCCQE_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFM_pCCQE_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFM_pCCQE_energy->Fill(aEvent.enu_truth,rweight);
						} else if ( aEvent.mode_truth == 1 ) {
							evt_tru0_recFM_pCCRES_count += ( 1. * rweight );
							evt_tru0_recFM_pCCRES_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFM_pCCRES_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFM_pCCRES_energy->Fill(aEvent.enu_truth,rweight);
						} else if ( aEvent.mode_truth == 2 ) {
							evt_tru0_recFM_pCCDIS_count += ( 1. * rweight );
							evt_tru0_recFM_pCCDIS_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFM_pCCDIS_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFM_pCCDIS_energy->Fill(aEvent.enu_truth,rweight);
						} else if ( aEvent.mode_truth == 3 ) {
							evt_tru0_recFM_pCCOTHER_count += ( 1. * rweight );
							evt_tru0_recFM_pCCOTHER_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFM_pCCOTHER_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFM_pCCOTHER_energy->Fill(aEvent.enu_truth,rweight);
						} // ccdis
						
					}

	
	
					if ( aEvent.trk_charge_minos < 0 ) {
	
						evt_tru0_recFMQ_p0_count += ( 1. * rweight );
						evt_tru0_recFMQ_p0_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
						evt_tru0_recFMQ_p0_momentum->Fill(aEvent.lep_mom_truth,rweight);
						evt_tru0_recFMQ_p0_energy->Fill(aEvent.enu_truth,rweight);
						if ( aEvent.ccnc_truth != 0 || aEvent.nuPDG_truth != 14 ) {
							evt_tru0_recFMQ_pBKG_count += ( 1. * rweight );
							evt_tru0_recFMQ_pBKG_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFMQ_pBKG_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFMQ_pBKG_energy->Fill(aEvent.enu_truth,rweight);									
							if ( aEvent.ccnc_truth == 1 ) { //NC
								evt_tru0_recFMQ_pNC_count += ( 1. * rweight );
								evt_tru0_recFMQ_pNC_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
								evt_tru0_recFMQ_pNC_momentum->Fill(aEvent.lep_mom_truth,rweight);
								evt_tru0_recFMQ_pNC_energy->Fill(aEvent.enu_truth,rweight);									
							}	else if ( aEvent.nuPDG_truth == -14 ) { // not NC and not a muon neutrino (probably anti neutrino??)
								evt_tru0_recFMQ_pWS_count += ( 1. * rweight ); 
								evt_tru0_recFMQ_pWS_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
								evt_tru0_recFMQ_pWS_momentum->Fill(aEvent.lep_mom_truth,rweight);
								evt_tru0_recFMQ_pWS_energy->Fill(aEvent.enu_truth,rweight);
							}	else if ( aEvent.nuPDG_truth != 14 ) { // not NC and not a muon neutrino (probably anti neutrino??)
								evt_tru0_recFMQ_pBOTHER_count += ( 1. * rweight ); 
								evt_tru0_recFMQ_pBOTHER_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
								evt_tru0_recFMQ_pBOTHER_momentum->Fill(aEvent.lep_mom_truth,rweight);
								evt_tru0_recFMQ_pBOTHER_energy->Fill(aEvent.enu_truth,rweight);
							}
						} else if ( aEvent.ccnc_truth == 0 && aEvent.nuPDG_truth == 14 ) {
		
							evt_tru0_recFMQ_pCC_count += ( 1. * rweight );
							evt_tru0_recFMQ_pCC_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
							evt_tru0_recFMQ_pCC_momentum->Fill(aEvent.lep_mom_truth,rweight);
							evt_tru0_recFMQ_pCC_energy->Fill(aEvent.enu_truth,rweight);
							if ( aEvent.lep_mom_truth < 25. ) 									evt_tru0_recFMQ_pCC_count_mom   += ( 1. * rweight );
							if ( (180./3.14159)*TMath::ACos(dcosztruth) < 36. ) evt_tru0_recFMQ_pCC_count_theta += ( 1. * rweight );

							// interaction type
							if ( aEvent.mode_truth == 0 ) {
								evt_tru0_recFMQ_pCCQE_count += ( 1. * rweight ); 
								evt_tru0_recFMQ_pCCQE_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
								evt_tru0_recFMQ_pCCQE_momentum->Fill(aEvent.lep_mom_truth,rweight);
								evt_tru0_recFMQ_pCCQE_energy->Fill(aEvent.enu_truth,rweight);
							} else if ( aEvent.mode_truth == 1 ) {
								evt_tru0_recFMQ_pCCRES_count += ( 1. * rweight );
								evt_tru0_recFMQ_pCCRES_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
								evt_tru0_recFMQ_pCCRES_momentum->Fill(aEvent.lep_mom_truth,rweight);
								evt_tru0_recFMQ_pCCRES_energy->Fill(aEvent.enu_truth,rweight);
							} else if ( aEvent.mode_truth == 2 ) {
								evt_tru0_recFMQ_pCCDIS_count += ( 1. * rweight );
								evt_tru0_recFMQ_pCCDIS_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
								evt_tru0_recFMQ_pCCDIS_momentum->Fill(aEvent.lep_mom_truth,rweight);
								evt_tru0_recFMQ_pCCDIS_energy->Fill(aEvent.enu_truth,rweight);
							} else if ( aEvent.mode_truth == 3 ) {
								evt_tru0_recFMQ_pCCOTHER_count += ( 1. * rweight );
								evt_tru0_recFMQ_pCCOTHER_theta->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
								evt_tru0_recFMQ_pCCOTHER_momentum->Fill(aEvent.lep_mom_truth,rweight);
								evt_tru0_recFMQ_pCCOTHER_energy->Fill(aEvent.enu_truth,rweight);
							}
						}
					}
				} else if ( rdiff < 1.1*RDIFF_CUT || thetadiff < 1.1*THETADIFF_CUT ) {
					// want to put some kind of measurement here about the best fit...
				}
			}
		}
	}

	cout << endl;
	cout << endl;
  cout << "                         Efficiency\t\tReconstructed\tPossible" << endl;
	cout << "-----------------------------------------------------------------------------------" << endl;
	cout << "      Total              "   << evt_tru0_recFMQ_pCC_count       / evt_truF_rec0_pCC_count       << "\t\t" << evt_tru0_recFMQ_pCC_count       << "\t\t" << evt_truF_rec0_pCC_count       << endl;
	cout << "         CCQE              " << evt_tru0_recFMQ_pCCQE_count     / evt_truF_rec0_pCCQE_count     << "\t\t" << evt_tru0_recFMQ_pCCQE_count     << "\t\t" << evt_truF_rec0_pCCQE_count     << endl;
	cout << "        CCRES              "	<< evt_tru0_recFMQ_pCCRES_count    / evt_truF_rec0_pCCRES_count    << "\t\t" << evt_tru0_recFMQ_pCCRES_count    << "\t\t" << evt_truF_rec0_pCCRES_count    << endl;
	cout << "        CCDIS              "	<< evt_tru0_recFMQ_pCCDIS_count    / evt_truF_rec0_pCCDIS_count    << "\t\t" << evt_tru0_recFMQ_pCCDIS_count    << "\t\t" << evt_truF_rec0_pCCDIS_count    << endl;
	cout << "      CCOTHER              "	<< evt_tru0_recFMQ_pCCOTHER_count  / evt_truF_rec0_pCCOTHER_count  << "\t\t" << evt_tru0_recFMQ_pCCOTHER_count  << "\t\t" << evt_truF_rec0_pCCOTHER_count  << endl;
	cout << "      Angle              "   << evt_tru0_recFMQ_pCC_count_theta / evt_truF_rec0_pCC_count_theta << "\t\t" << evt_tru0_recFMQ_pCC_count_theta << "\t\t" << evt_truF_rec0_pCC_count_theta << endl;
	cout << "   Momentum              "   << evt_tru0_recFMQ_pCC_count_mom   / evt_truF_rec0_pCC_count_mom   << "\t\t" << evt_tru0_recFMQ_pCC_count_mom   << "\t\t" << evt_truF_rec0_pCC_count_mom   << endl;
	cout << "-----------------------------------------------------------------------------------" << endl;
	cout << "   Fiducial (reco)       "   << evt_tru0_recF_p0_count       / evt_truF_rec0_pCC_count      << "\t\t" << evt_tru0_recF_p0_count       << "\t\t" << evt_truF_rec0_pCC_count      << endl;
	cout << "         CCQE              " << evt_tru0_recF_pCCQE_count    / evt_truF_rec0_pCCQE_count    << "\t\t" << evt_tru0_recF_pCCQE_count    << "\t\t" << evt_truF_rec0_pCCQE_count    << endl;
	cout << "        CCRES              " << evt_tru0_recF_pCCRES_count   / evt_truF_rec0_pCCRES_count   << "\t\t" << evt_tru0_recF_pCCRES_count   << "\t\t" << evt_truF_rec0_pCCRES_count   << endl;
	cout << "        CCDIS              " << evt_tru0_recF_pCCDIS_count   / evt_truF_rec0_pCCDIS_count   << "\t\t" << evt_tru0_recF_pCCDIS_count   << "\t\t" << evt_truF_rec0_pCCDIS_count   << endl;
	cout << "      CCOTHER              " << evt_tru0_recF_pCCOTHER_count / evt_truF_rec0_pCCOTHER_count << "\t\t" << evt_tru0_recF_pCCOTHER_count << "\t\t" << evt_truF_rec0_pCCOTHER_count << endl;
	cout << "-----------------------------------------------------------------------------------" << endl;
	cout << "      Match (reco)       "   << evt_tru0_recFM_p0_count       / evt_tru0_recF_pCC_count      << "\t\t" << evt_tru0_recFM_p0_count       << "\t\t" << evt_tru0_recF_pCC_count      << endl;
	cout << "         CCQE              " << evt_tru0_recFM_pCCQE_count    / evt_tru0_recF_pCCQE_count    << "\t\t" << evt_tru0_recFM_pCCQE_count    << "\t\t" << evt_tru0_recF_pCCQE_count    << endl;
	cout << "        CCRES              " << evt_tru0_recFM_pCCRES_count   / evt_tru0_recF_pCCRES_count   << "\t\t" << evt_tru0_recFM_pCCRES_count   << "\t\t" << evt_tru0_recF_pCCRES_count   << endl;
	cout << "        CCDIS              " << evt_tru0_recFM_pCCDIS_count   / evt_tru0_recF_pCCDIS_count   << "\t\t" << evt_tru0_recFM_pCCDIS_count   << "\t\t" << evt_tru0_recF_pCCDIS_count   << endl;
	cout << "      CCOTHER              " << evt_tru0_recFM_pCCOTHER_count / evt_tru0_recF_pCCOTHER_count << "\t\t" << evt_tru0_recFM_pCCOTHER_count << "\t\t" << evt_tru0_recF_pCCOTHER_count << endl;
	cout << "-----------------------------------------------------------------------------------" << endl;
	cout << "     Charge (reco)       "   << evt_tru0_recFMQ_p0_count       / evt_tru0_recF_pCC_count      << "\t\t" << evt_tru0_recFMQ_p0_count       << "\t\t" << evt_tru0_recF_pCC_count      << endl;
	cout << "         CCQE              " << evt_tru0_recFMQ_pCCQE_count    / evt_tru0_recF_pCCQE_count    << "\t\t" << evt_tru0_recFMQ_pCCQE_count    << "\t\t" << evt_tru0_recF_pCCQE_count    << endl;
	cout << "        CCRES              " << evt_tru0_recFMQ_pCCRES_count   / evt_tru0_recF_pCCRES_count   << "\t\t" << evt_tru0_recFMQ_pCCRES_count   << "\t\t" << evt_tru0_recF_pCCRES_count   << endl;
	cout << "        CCDIS              " << evt_tru0_recFMQ_pCCDIS_count   / evt_tru0_recF_pCCDIS_count   << "\t\t" << evt_tru0_recFMQ_pCCDIS_count   << "\t\t" << evt_tru0_recF_pCCDIS_count   << endl;
	cout << "      CCOTHER              " << evt_tru0_recFMQ_pCCOTHER_count / evt_tru0_recF_pCCOTHER_count << "\t\t" << evt_tru0_recFMQ_pCCOTHER_count << "\t\t" << evt_tru0_recF_pCCOTHER_count << endl;
	cout << "-----------------------------------------------------------------------------------" << endl;
	cout << "      MINOS (true)       "   << evt_truFQ_rec0_pCC_count      / evt_truF_rec0_pCC_count      << "\t\t" << evt_truFQ_rec0_pCC_count      << "\t\t" << evt_truF_rec0_pCC_count      << endl;
	cout << "         CCQE              " << evt_truFQ_rec0_pCCQE_count    / evt_truF_rec0_pCCQE_count    << "\t\t" << evt_truFQ_rec0_pCCQE_count    << "\t\t" << evt_truF_rec0_pCCQE_count    << endl;
	cout << "        CCRES              " << evt_truFQ_rec0_pCCRES_count   / evt_truF_rec0_pCCRES_count   << "\t\t" << evt_truFQ_rec0_pCCRES_count   << "\t\t" << evt_truF_rec0_pCCRES_count   << endl;
	cout << "        CCDIS              " << evt_truFQ_rec0_pCCDIS_count   / evt_truF_rec0_pCCDIS_count   << "\t\t" << evt_truFQ_rec0_pCCDIS_count   << "\t\t" << evt_truF_rec0_pCCDIS_count   << endl;
	cout << "      CCOTHER              " << evt_truFQ_rec0_pCCOTHER_count / evt_truF_rec0_pCCOTHER_count << "\t\t" << evt_truFQ_rec0_pCCOTHER_count << "\t\t" << evt_truF_rec0_pCCOTHER_count << endl;
  cout << "-----------------------------------------------------------------------------------" << endl;
  cout << endl;
  cout << endl;
  cout << "                          Purity\t\tCount\t\tTotal" << endl;
  cout << "-----------------------------------------------------------------------------------"           << endl;
  cout << "     Fiducial            "   <<   evt_tru0_recF_pCC_count      /   evt_tru0_recF_p0_count << "\t\t" << evt_tru0_recF_pCC_count      << "\t\t" << evt_tru0_recF_p0_count << endl;
  cout << "           CCQE            " <<   evt_tru0_recF_pCCQE_count    /   evt_tru0_recF_p0_count << "\t\t" << evt_tru0_recF_pCCQE_count    << endl;
  cout << "          CCDIS            " <<   evt_tru0_recF_pCCDIS_count   /   evt_tru0_recF_p0_count << "\t\t" << evt_tru0_recF_pCCDIS_count   << endl;
  cout << "          CCRES            " <<   evt_tru0_recF_pCCRES_count   /   evt_tru0_recF_p0_count << "\t\t" << evt_tru0_recF_pCCRES_count   << endl;
  cout << "        CCOTHER            " <<   evt_tru0_recF_pCCOTHER_count /   evt_tru0_recF_p0_count << "\t\t" << evt_tru0_recF_pCCOTHER_count << endl;
  cout << "             NC            " <<   evt_tru0_recF_pNC_count      /   evt_tru0_recF_p0_count << "\t\t" << evt_tru0_recF_pNC_count      << endl;
  cout << "             WS            " <<   evt_tru0_recF_pWS_count      /   evt_tru0_recF_p0_count << "\t\t" << evt_tru0_recF_pWS_count      << endl;
  cout << "          OTHER            " <<   evt_tru0_recF_pBOTHER_count  /   evt_tru0_recF_p0_count << "\t\t" << evt_tru0_recF_pBOTHER_count  << endl;
  cout << "-----------------------------------------------------------------------------------"           << endl;
  cout << "        Match            "   <<   evt_tru0_recFM_pCC_count      /   evt_tru0_recFM_p0_count << "\t\t" << evt_tru0_recFM_pCC_count      << "\t\t" << evt_tru0_recFM_p0_count << endl;
  cout << "           CCQE            " <<   evt_tru0_recFM_pCCQE_count    /   evt_tru0_recFM_p0_count << "\t\t" << evt_tru0_recFM_pCCQE_count    << endl;
  cout << "          CCDIS            " <<   evt_tru0_recFM_pCCDIS_count   /   evt_tru0_recFM_p0_count << "\t\t" << evt_tru0_recFM_pCCDIS_count   << endl;
  cout << "          CCRES            " <<   evt_tru0_recFM_pCCRES_count   /   evt_tru0_recFM_p0_count << "\t\t" << evt_tru0_recFM_pCCRES_count   << endl;
  cout << "        CCOTHER            " <<   evt_tru0_recFM_pCCOTHER_count /   evt_tru0_recFM_p0_count << "\t\t" << evt_tru0_recFM_pCCOTHER_count << endl;
  cout << "             NC            " <<   evt_tru0_recFM_pNC_count      /   evt_tru0_recFM_p0_count << "\t\t" << evt_tru0_recFM_pNC_count      << endl;
  cout << "             WS            " <<   evt_tru0_recFM_pWS_count      /   evt_tru0_recFM_p0_count << "\t\t" << evt_tru0_recFM_pWS_count      << endl;
  cout << "          OTHER            " <<   evt_tru0_recFM_pBOTHER_count  /   evt_tru0_recFM_p0_count << "\t\t" << evt_tru0_recFM_pBOTHER_count  << endl;
  cout << "-----------------------------------------------------------------------------------" << endl;
  cout << "       Charge             "  <<   evt_tru0_recFMQ_pCC_count      /   evt_tru0_recFMQ_p0_count << "\t\t" << evt_tru0_recFMQ_pCC_count      << "\t\t" << evt_tru0_recFMQ_p0_count << endl;
  cout << "           CCQE            " <<   evt_tru0_recFMQ_pCCQE_count    /   evt_tru0_recFMQ_p0_count << "\t\t" << evt_tru0_recFMQ_pCCQE_count    << endl;
  cout << "          CCDIS            " <<   evt_tru0_recFMQ_pCCDIS_count   /   evt_tru0_recFMQ_p0_count << "\t\t" << evt_tru0_recFMQ_pCCDIS_count   << endl;
  cout << "          CCRES            " <<   evt_tru0_recFMQ_pCCRES_count   /   evt_tru0_recFMQ_p0_count << "\t\t" << evt_tru0_recFMQ_pCCRES_count   << endl;
  cout << "        CCOTHER            " <<   evt_tru0_recFMQ_pCCOTHER_count /   evt_tru0_recFMQ_p0_count << "\t\t" << evt_tru0_recFMQ_pCCOTHER_count << endl;
  cout << "             NC            " <<   evt_tru0_recFMQ_pNC_count      /   evt_tru0_recFMQ_p0_count << "\t\t" << evt_tru0_recFMQ_pNC_count      << endl;
  cout << "             WS            " <<   evt_tru0_recFMQ_pWS_count      /   evt_tru0_recFMQ_p0_count << "\t\t" << evt_tru0_recFMQ_pWS_count      << endl;
  cout << "          OTHER            " <<   evt_tru0_recFMQ_pBOTHER_count  /   evt_tru0_recFMQ_p0_count << "\t\t" << evt_tru0_recFMQ_pBOTHER_count  << endl;
  cout << "-----------------------------------------------------------------------------------"           << endl;
  
  
  cout << endl;
  cout << endl;
  cout << "Drawing plots..." << endl;


	new TCanvas;	
	{
		TH1F *count_truF_pNU_energy = evt_truF_rec0_pNU_energy->Clone("count_truF_pNU_energy");
		TH1F *count_truF_pANU_energy = evt_truF_rec0_pANU_energy->Clone("count_truF_pANU_energy");
    count_truF_pNU_energy->SetTitle("True energy in true fiducial volume");
		count_truF_pNU_energy->SetXTitle("E_{#mu} (GeV)");
		count_truF_pNU_energy->SetLineColor(33);
		count_truF_pNU_energy->Draw();
		count_truF_pANU_energy->SetLineColor(31);
		count_truF_pANU_energy->Draw("same");
 		count_truF_pNUANU_energy_leg = new TLegend(0.6,0.2,0.79,0.39);
		count_truF_pNUANU_energy_leg->AddEntry(count_truF_pNU_energy,"#nu","l");
		count_truF_pNUANU_energy_leg->AddEntry(count_truF_pANU_energy,"anti-#nu","l");
		count_truF_pNUANU_energy_leg->Draw("same");
	} gPad->GetCanvas()->Print("count_truF_pNUANU_energy.png");

	// efficiencies
	new TCanvas;	
	{
		TH1F *efficiency_recFMQ_pCC_theta = evt_tru0_recFMQ_pCC_theta->Clone("efficiency_recFMQ_pCC_theta");
			efficiency_recFMQ_pCC_theta->Divide(evt_truF_rec0_pCC_theta);
      efficiency_recFMQ_pCC_theta->SetTitle("Efficiency of all cuts (recFMQ/truF)");
			efficiency_recFMQ_pCC_theta->SetXTitle("#theta_{#mu} (degrees)");
			efficiency_recFMQ_pCC_theta->SetYTitle("efficiency");
			efficiency_recFMQ_pCC_theta->SetMinimum(0);
			efficiency_recFMQ_pCC_theta->SetMaximum(1);
			efficiency_recFMQ_pCC_theta->SetLineColor(31);
			efficiency_recFMQ_pCC_theta->SetStats(0);
		efficiency_recFMQ_pCC_theta->Draw();
	} gPad->GetCanvas()->Print("efficiency_recFMQ_pCC_theta.png");
	{
		TH1F *efficiency_recFMQ_pCC_momentum = evt_tru0_recFMQ_pCC_momentum->Clone("efficiency_recFMQ_pCC_momentum");
			efficiency_recFMQ_pCC_momentum->Divide(evt_truF_rec0_pCC_momentum);
      efficiency_recFMQ_pCC_momentum->SetTitle("Efficiency of all cuts (recFMQ/truF)");
			efficiency_recFMQ_pCC_momentum->SetXTitle("p_{#mu} (GeV/c)");
			efficiency_recFMQ_pCC_momentum->SetYTitle("efficiency");
			efficiency_recFMQ_pCC_momentum->SetMinimum(0);
			efficiency_recFMQ_pCC_momentum->SetMaximum(1);
			efficiency_recFMQ_pCC_momentum->SetLineColor(31);
			efficiency_recFMQ_pCC_momentum->SetStats(0);
		efficiency_recFMQ_pCC_momentum->Draw();
	} gPad->GetCanvas()->Print("efficiency_recFMQ_pCC_momentum.png");

	// purities
	{
	int purity_pCC_color       = 31;
	int purity_pCCQE_color     = 30;
	int purity_pCCRES_color    = 31;
	int purity_pCCDIS_color    = 32;
	int purity_pCCOTHER_color  = 33;
	int purity_pBKG_color      = 20;
	int purity_pNC_color       = 20;
	int purity_pWS_color       = 24;
	int purity_pBOTHER_color   = 27;

	TH1F *purity_tru0_recF_pCC_theta = evt_tru0_recF_pCC_theta->Clone("purity_tru0_recF_pCC_theta");
	TH1F *purity_tru0_recF_pCC_momentum = evt_tru0_recF_pCC_momentum->Clone("purity_tru0_recF_pCC_momentum");
	TH1F *purity_tru0_recFM_pCC_theta = evt_tru0_recFM_pCC_theta->Clone("purity_tru0_recFM_pCC_theta");
	TH1F *purity_tru0_recFM_pCC_momentum = evt_tru0_recFM_pCC_momentum->Clone("purity_tru0_recFM_pCC_momentum");
	TH1F *purity_tru0_recFMQ_pCC_theta = evt_tru0_recFMQ_pCC_theta->Clone("purity_tru0_recFMQ_pCC_theta");
	TH1F *purity_tru0_recFMQ_pCC_momentum = evt_tru0_recFMQ_pCC_momentum->Clone("purity_tru0_recFMQ_pCC_momentum");
		purity_tru0_recF_pCC_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recF_pCC_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recFM_pCC_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFM_pCC_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFMQ_pCC_theta->Divide(evt_tru0_recFMQ_p0_theta);
		purity_tru0_recFMQ_pCC_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
			purity_tru0_recF_pCC_theta->SetLineWidth(1);
			purity_tru0_recF_pCC_theta->SetLineColor(purity_pCC_color);
			purity_tru0_recF_pCC_theta->SetFillColor(purity_pCC_color);
			purity_tru0_recF_pCC_momentum->SetLineWidth(1);
			purity_tru0_recF_pCC_momentum->SetLineColor(purity_pCC_color);
			purity_tru0_recF_pCC_momentum->SetFillColor(purity_pCC_color);
			purity_tru0_recFM_pCC_theta->SetLineWidth(1);
			purity_tru0_recFM_pCC_theta->SetLineColor(purity_pCC_color);
			purity_tru0_recFM_pCC_theta->SetFillColor(purity_pCC_color);
			purity_tru0_recFM_pCC_momentum->SetLineWidth(1);
			purity_tru0_recFM_pCC_momentum->SetLineColor(purity_pCC_color);
			purity_tru0_recFM_pCC_momentum->SetFillColor(purity_pCC_color);
			purity_tru0_recFMQ_pCC_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pCC_theta->SetLineColor(purity_pCC_color);
			purity_tru0_recFMQ_pCC_theta->SetFillColor(purity_pCC_color);
			purity_tru0_recFMQ_pCC_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pCC_momentum->SetLineColor(purity_pCC_color);
			purity_tru0_recFMQ_pCC_momentum->SetFillColor(purity_pCC_color);

	TH1F *purity_tru0_recF_pCCDIS_theta = evt_tru0_recF_pCCDIS_theta->Clone("purity_tru0_recF_pCCDIS_theta");
	TH1F *purity_tru0_recF_pCCDIS_momentum = evt_tru0_recF_pCCDIS_momentum->Clone("purity_tru0_recF_pCCDIS_momentum");
	TH1F *purity_tru0_recFM_pCCDIS_theta = evt_tru0_recFM_pCCDIS_theta->Clone("purity_tru0_recFM_pCCDIS_theta");
	TH1F *purity_tru0_recFM_pCCDIS_momentum = evt_tru0_recFM_pCCDIS_momentum->Clone("purity_tru0_recFM_pCCDIS_momentum");
	TH1F *purity_tru0_recFMQ_pCCDIS_theta = evt_tru0_recFMQ_pCCDIS_theta->Clone("purity_tru0_recFMQ_pCCDIS_theta");
	TH1F *purity_tru0_recFMQ_pCCDIS_momentum = evt_tru0_recFMQ_pCCDIS_momentum->Clone("purity_tru0_recFMQ_pCCDIS_momentum");
		purity_tru0_recF_pCCDIS_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recF_pCCDIS_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recFM_pCCDIS_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFM_pCCDIS_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFMQ_pCCDIS_theta->Divide(evt_tru0_recFMQ_p0_theta);
		purity_tru0_recFMQ_pCCDIS_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
			purity_tru0_recF_pCCDIS_theta->SetLineWidth(1);
			purity_tru0_recF_pCCDIS_theta->SetLineColor(purity_pCCDIS_color);
			purity_tru0_recF_pCCDIS_theta->SetFillColor(purity_pCCDIS_color);
			purity_tru0_recF_pCCDIS_momentum->SetLineWidth(1);
			purity_tru0_recF_pCCDIS_momentum->SetLineColor(purity_pCCDIS_color);
			purity_tru0_recF_pCCDIS_momentum->SetFillColor(purity_pCCDIS_color);
			purity_tru0_recFM_pCCDIS_theta->SetLineWidth(1);
			purity_tru0_recFM_pCCDIS_theta->SetLineColor(purity_pCCDIS_color);
			purity_tru0_recFM_pCCDIS_theta->SetFillColor(purity_pCCDIS_color);
			purity_tru0_recFM_pCCDIS_momentum->SetLineWidth(1);
			purity_tru0_recFM_pCCDIS_momentum->SetLineColor(purity_pCCDIS_color);
			purity_tru0_recFM_pCCDIS_momentum->SetFillColor(purity_pCCDIS_color);
			purity_tru0_recFMQ_pCCDIS_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pCCDIS_theta->SetLineColor(purity_pCCDIS_color);
			purity_tru0_recFMQ_pCCDIS_theta->SetFillColor(purity_pCCDIS_color);
			purity_tru0_recFMQ_pCCDIS_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pCCDIS_momentum->SetLineColor(purity_pCCDIS_color);
			purity_tru0_recFMQ_pCCDIS_momentum->SetFillColor(purity_pCCDIS_color);

	TH1F *purity_tru0_recF_pCCRES_theta = evt_tru0_recF_pCCRES_theta->Clone("purity_tru0_recF_pCCRES_theta");
	TH1F *purity_tru0_recF_pCCRES_momentum = evt_tru0_recF_pCCRES_momentum->Clone("purity_tru0_recF_pCCRES_momentum");
	TH1F *purity_tru0_recFM_pCCRES_theta = evt_tru0_recFM_pCCRES_theta->Clone("purity_tru0_recFM_pCCRES_theta");
	TH1F *purity_tru0_recFM_pCCRES_momentum = evt_tru0_recFM_pCCRES_momentum->Clone("purity_tru0_recFM_pCCRES_momentum");
	TH1F *purity_tru0_recFMQ_pCCRES_theta = evt_tru0_recFMQ_pCCRES_theta->Clone("purity_tru0_recFMQ_pCCRES_theta");
	TH1F *purity_tru0_recFMQ_pCCRES_momentum = evt_tru0_recFMQ_pCCRES_momentum->Clone("purity_tru0_recFMQ_pCCRES_momentum");
		purity_tru0_recF_pCCRES_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recF_pCCRES_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recFM_pCCRES_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFM_pCCRES_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFMQ_pCCRES_theta->Divide(evt_tru0_recFMQ_p0_theta);
		purity_tru0_recFMQ_pCCRES_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
			purity_tru0_recF_pCCRES_theta->SetLineWidth(1);
			purity_tru0_recF_pCCRES_theta->SetLineColor(purity_pCCRES_color);
			purity_tru0_recF_pCCRES_theta->SetFillColor(purity_pCCRES_color);
			purity_tru0_recF_pCCRES_momentum->SetLineWidth(1);
			purity_tru0_recF_pCCRES_momentum->SetLineColor(purity_pCCRES_color);
			purity_tru0_recF_pCCRES_momentum->SetFillColor(purity_pCCRES_color);
			purity_tru0_recFM_pCCRES_theta->SetLineWidth(1);
			purity_tru0_recFM_pCCRES_theta->SetLineColor(purity_pCCRES_color);
			purity_tru0_recFM_pCCRES_theta->SetFillColor(purity_pCCRES_color);
			purity_tru0_recFM_pCCRES_momentum->SetLineWidth(1);
			purity_tru0_recFM_pCCRES_momentum->SetLineColor(purity_pCCRES_color);
			purity_tru0_recFM_pCCRES_momentum->SetFillColor(purity_pCCRES_color);
			purity_tru0_recFMQ_pCCRES_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pCCRES_theta->SetLineColor(purity_pCCRES_color);
			purity_tru0_recFMQ_pCCRES_theta->SetFillColor(purity_pCCRES_color);
			purity_tru0_recFMQ_pCCRES_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pCCRES_momentum->SetLineColor(purity_pCCRES_color);
			purity_tru0_recFMQ_pCCRES_momentum->SetFillColor(purity_pCCRES_color);

	TH1F *purity_tru0_recF_pCCQE_theta = evt_tru0_recF_pCCQE_theta->Clone("purity_tru0_recF_pCCQE_theta");
	TH1F *purity_tru0_recF_pCCQE_momentum = evt_tru0_recF_pCCQE_momentum->Clone("purity_tru0_recF_pCCQE_momentum");
	TH1F *purity_tru0_recFM_pCCQE_theta = evt_tru0_recFM_pCCQE_theta->Clone("purity_tru0_recFM_pCCQE_theta");
	TH1F *purity_tru0_recFM_pCCQE_momentum = evt_tru0_recFM_pCCQE_momentum->Clone("purity_tru0_recFM_pCCQE_momentum");
	TH1F *purity_tru0_recFMQ_pCCQE_theta = evt_tru0_recFMQ_pCCQE_theta->Clone("purity_tru0_recFMQ_pCCQE_theta");
	TH1F *purity_tru0_recFMQ_pCCQE_momentum = evt_tru0_recFMQ_pCCQE_momentum->Clone("purity_tru0_recFMQ_pCCQE_momentum");
		purity_tru0_recF_pCCQE_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recF_pCCQE_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recFM_pCCQE_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFM_pCCQE_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFMQ_pCCQE_theta->Divide(evt_tru0_recFMQ_p0_theta);
		purity_tru0_recFMQ_pCCQE_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
			purity_tru0_recF_pCCQE_theta->SetLineWidth(1);
			purity_tru0_recF_pCCQE_theta->SetLineColor(purity_pCCQE_color);
			purity_tru0_recF_pCCQE_theta->SetFillColor(purity_pCCQE_color);
			purity_tru0_recF_pCCQE_momentum->SetLineWidth(1);
			purity_tru0_recF_pCCQE_momentum->SetLineColor(purity_pCCQE_color);
			purity_tru0_recF_pCCQE_momentum->SetFillColor(purity_pCCQE_color);
			purity_tru0_recFM_pCCQE_theta->SetLineWidth(1);
			purity_tru0_recFM_pCCQE_theta->SetLineColor(purity_pCCQE_color);
			purity_tru0_recFM_pCCQE_theta->SetFillColor(purity_pCCQE_color);
			purity_tru0_recFM_pCCQE_momentum->SetLineWidth(1);
			purity_tru0_recFM_pCCQE_momentum->SetLineColor(purity_pCCQE_color);
			purity_tru0_recFM_pCCQE_momentum->SetFillColor(purity_pCCQE_color);
			purity_tru0_recFMQ_pCCQE_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pCCQE_theta->SetLineColor(purity_pCCQE_color);
			purity_tru0_recFMQ_pCCQE_theta->SetFillColor(purity_pCCQE_color);
			purity_tru0_recFMQ_pCCQE_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pCCQE_momentum->SetLineColor(purity_pCCQE_color);
			purity_tru0_recFMQ_pCCQE_momentum->SetFillColor(purity_pCCQE_color);

	TH1F *purity_tru0_recF_pCCOTHER_theta = evt_tru0_recF_pCCOTHER_theta->Clone("purity_tru0_recF_pCCOTHER_theta");
	TH1F *purity_tru0_recF_pCCOTHER_momentum = evt_tru0_recF_pCCOTHER_momentum->Clone("purity_tru0_recF_pCCOTHER_momentum");
	TH1F *purity_tru0_recFM_pCCOTHER_theta = evt_tru0_recFM_pCCOTHER_theta->Clone("purity_tru0_recFM_pCCOTHER_theta");
	TH1F *purity_tru0_recFM_pCCOTHER_momentum = evt_tru0_recFM_pCCOTHER_momentum->Clone("purity_tru0_recFM_pCCOTHER_momentum");
	TH1F *purity_tru0_recFMQ_pCCOTHER_theta = evt_tru0_recFMQ_pCCOTHER_theta->Clone("purity_tru0_recFMQ_pCCOTHER_theta");
	TH1F *purity_tru0_recFMQ_pCCOTHER_momentum = evt_tru0_recFMQ_pCCOTHER_momentum->Clone("purity_tru0_recFMQ_pCCOTHER_momentum");
		purity_tru0_recF_pCCOTHER_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recF_pCCOTHER_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recFM_pCCOTHER_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFM_pCCOTHER_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFMQ_pCCOTHER_theta->Divide(evt_tru0_recFMQ_p0_theta);
		purity_tru0_recFMQ_pCCOTHER_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
			purity_tru0_recF_pCCOTHER_theta->SetLineWidth(1);
			purity_tru0_recF_pCCOTHER_theta->SetLineColor(purity_pCCOTHER_color);
			purity_tru0_recF_pCCOTHER_theta->SetFillColor(purity_pCCOTHER_color);
			purity_tru0_recF_pCCOTHER_momentum->SetLineWidth(1);
			purity_tru0_recF_pCCOTHER_momentum->SetLineColor(purity_pCCOTHER_color);
			purity_tru0_recF_pCCOTHER_momentum->SetFillColor(purity_pCCOTHER_color);
			purity_tru0_recFM_pCCOTHER_theta->SetLineWidth(1);
			purity_tru0_recFM_pCCOTHER_theta->SetLineColor(purity_pCCOTHER_color);
			purity_tru0_recFM_pCCOTHER_theta->SetFillColor(purity_pCCOTHER_color);
			purity_tru0_recFM_pCCOTHER_momentum->SetLineWidth(1);
			purity_tru0_recFM_pCCOTHER_momentum->SetLineColor(purity_pCCOTHER_color);
			purity_tru0_recFM_pCCOTHER_momentum->SetFillColor(purity_pCCOTHER_color);
			purity_tru0_recFMQ_pCCOTHER_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pCCOTHER_theta->SetLineColor(purity_pCCOTHER_color);
			purity_tru0_recFMQ_pCCOTHER_theta->SetFillColor(purity_pCCOTHER_color);
			purity_tru0_recFMQ_pCCOTHER_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pCCOTHER_momentum->SetLineColor(purity_pCCOTHER_color);
			purity_tru0_recFMQ_pCCOTHER_momentum->SetFillColor(purity_pCCOTHER_color);

	TH1F *purity_tru0_recF_pBKG_theta = evt_tru0_recF_pBKG_theta->Clone("purity_tru0_recF_pBKG_theta");
	TH1F *purity_tru0_recF_pBKG_momentum = evt_tru0_recF_pBKG_momentum->Clone("purity_tru0_recF_pBKG_momentum");
	TH1F *purity_tru0_recFM_pBKG_theta = evt_tru0_recFM_pBKG_theta->Clone("purity_tru0_recFM_pBKG_theta");
	TH1F *purity_tru0_recFM_pBKG_momentum = evt_tru0_recFM_pBKG_momentum->Clone("purity_tru0_recFM_pBKG_momentum");
	TH1F *purity_tru0_recFMQ_pBKG_theta = evt_tru0_recFMQ_pBKG_theta->Clone("purity_tru0_recFMQ_pBKG_theta");
	TH1F *purity_tru0_recFMQ_pBKG_momentum = evt_tru0_recFMQ_pBKG_momentum->Clone("purity_tru0_recFMQ_pBKG_momentum");
		purity_tru0_recF_pBKG_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recF_pBKG_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recFM_pBKG_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFM_pBKG_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFMQ_pBKG_theta->Divide(evt_tru0_recFMQ_p0_theta);
		purity_tru0_recFMQ_pBKG_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
			purity_tru0_recF_pBKG_theta->SetLineWidth(1);
			purity_tru0_recF_pBKG_theta->SetLineColor(purity_pBKG_color);
			purity_tru0_recF_pBKG_theta->SetFillColor(purity_pBKG_color);
			purity_tru0_recF_pBKG_momentum->SetLineWidth(1);
			purity_tru0_recF_pBKG_momentum->SetLineColor(purity_pBKG_color);
			purity_tru0_recF_pBKG_momentum->SetFillColor(purity_pBKG_color);
			purity_tru0_recFM_pBKG_theta->SetLineWidth(1);
			purity_tru0_recFM_pBKG_theta->SetLineColor(purity_pBKG_color);
			purity_tru0_recFM_pBKG_theta->SetFillColor(purity_pBKG_color);
			purity_tru0_recFM_pBKG_momentum->SetLineWidth(1);
			purity_tru0_recFM_pBKG_momentum->SetLineColor(purity_pBKG_color);
			purity_tru0_recFM_pBKG_momentum->SetFillColor(purity_pBKG_color);
			purity_tru0_recFMQ_pBKG_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pBKG_theta->SetLineColor(purity_pBKG_color);
			purity_tru0_recFMQ_pBKG_theta->SetFillColor(purity_pBKG_color);
			purity_tru0_recFMQ_pBKG_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pBKG_momentum->SetLineColor(purity_pBKG_color);
			purity_tru0_recFMQ_pBKG_momentum->SetFillColor(purity_pBKG_color);

	TH1F *purity_tru0_recF_pNC_theta = evt_tru0_recF_pNC_theta->Clone("purity_tru0_recF_pNC_theta");
	TH1F *purity_tru0_recF_pNC_momentum = evt_tru0_recF_pNC_momentum->Clone("purity_tru0_recF_pNC_momentum");
	TH1F *purity_tru0_recFM_pNC_theta = evt_tru0_recFM_pNC_theta->Clone("purity_tru0_recFM_pNC_theta");
	TH1F *purity_tru0_recFM_pNC_momentum = evt_tru0_recFM_pNC_momentum->Clone("purity_tru0_recFM_pNC_momentum");
	TH1F *purity_tru0_recFMQ_pNC_theta = evt_tru0_recFMQ_pNC_theta->Clone("purity_tru0_recFMQ_pNC_theta");
	TH1F *purity_tru0_recFMQ_pNC_momentum = evt_tru0_recFMQ_pNC_momentum->Clone("purity_tru0_recFMQ_pNC_momentum");
		purity_tru0_recF_pNC_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recF_pNC_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recFM_pNC_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFM_pNC_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFMQ_pNC_theta->Divide(evt_tru0_recFMQ_p0_theta);
		purity_tru0_recFMQ_pNC_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
			purity_tru0_recF_pNC_theta->SetLineWidth(1);
			purity_tru0_recF_pNC_theta->SetLineColor(purity_pNC_color);
			purity_tru0_recF_pNC_theta->SetFillColor(purity_pNC_color);
			purity_tru0_recF_pNC_momentum->SetLineWidth(1);
			purity_tru0_recF_pNC_momentum->SetLineColor(purity_pNC_color);
			purity_tru0_recF_pNC_momentum->SetFillColor(purity_pNC_color);
			purity_tru0_recFM_pNC_theta->SetLineWidth(1);
			purity_tru0_recFM_pNC_theta->SetLineColor(purity_pNC_color);
			purity_tru0_recFM_pNC_theta->SetFillColor(purity_pNC_color);
			purity_tru0_recFM_pNC_momentum->SetLineWidth(1);
			purity_tru0_recFM_pNC_momentum->SetLineColor(purity_pNC_color);
			purity_tru0_recFM_pNC_momentum->SetFillColor(purity_pNC_color);
			purity_tru0_recFMQ_pNC_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pNC_theta->SetLineColor(purity_pNC_color);
			purity_tru0_recFMQ_pNC_theta->SetFillColor(purity_pNC_color);
			purity_tru0_recFMQ_pNC_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pNC_momentum->SetLineColor(purity_pNC_color);
			purity_tru0_recFMQ_pNC_momentum->SetFillColor(purity_pNC_color);

	TH1F *purity_tru0_recF_pWS_theta = evt_tru0_recF_pWS_theta->Clone("purity_tru0_recF_pWS_theta");
	TH1F *purity_tru0_recF_pWS_momentum = evt_tru0_recF_pWS_momentum->Clone("purity_tru0_recF_pWS_momentum");
	TH1F *purity_tru0_recFM_pWS_theta = evt_tru0_recFM_pWS_theta->Clone("purity_tru0_recFM_pWS_theta");
	TH1F *purity_tru0_recFM_pWS_momentum = evt_tru0_recFM_pWS_momentum->Clone("purity_tru0_recFM_pWS_momentum");
	TH1F *purity_tru0_recFMQ_pWS_theta = evt_tru0_recFMQ_pWS_theta->Clone("purity_tru0_recFMQ_pWS_theta");
	TH1F *purity_tru0_recFMQ_pWS_momentum = evt_tru0_recFMQ_pWS_momentum->Clone("purity_tru0_recFMQ_pWS_momentum");
		purity_tru0_recF_pWS_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recF_pWS_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recFM_pWS_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFM_pWS_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFMQ_pWS_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
		purity_tru0_recFMQ_pWS_theta->Divide(evt_tru0_recFMQ_p0_theta);
			purity_tru0_recF_pWS_theta->SetLineWidth(1);
			purity_tru0_recF_pWS_theta->SetLineColor(purity_pWS_color);
			purity_tru0_recF_pWS_theta->SetFillColor(purity_pWS_color);
			purity_tru0_recF_pWS_momentum->SetLineWidth(1);
			purity_tru0_recF_pWS_momentum->SetLineColor(purity_pWS_color);
			purity_tru0_recF_pWS_momentum->SetFillColor(purity_pWS_color);
			purity_tru0_recFM_pWS_theta->SetLineWidth(1);
			purity_tru0_recFM_pWS_theta->SetLineColor(purity_pWS_color);
			purity_tru0_recFM_pWS_theta->SetFillColor(purity_pWS_color);
			purity_tru0_recFM_pWS_momentum->SetLineWidth(1);
			purity_tru0_recFM_pWS_momentum->SetLineColor(purity_pWS_color);
			purity_tru0_recFM_pWS_momentum->SetFillColor(purity_pWS_color);
			purity_tru0_recFMQ_pWS_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pWS_theta->SetLineColor(purity_pWS_color);
			purity_tru0_recFMQ_pWS_theta->SetFillColor(purity_pWS_color);
			purity_tru0_recFMQ_pWS_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pWS_momentum->SetLineColor(purity_pWS_color);
			purity_tru0_recFMQ_pWS_momentum->SetFillColor(purity_pWS_color);

	TH1F *purity_tru0_recF_pBOTHER_theta = evt_tru0_recF_pBOTHER_theta->Clone("purity_tru0_recF_pBOTHER_theta");
	TH1F *purity_tru0_recF_pBOTHER_momentum = evt_tru0_recF_pBOTHER_momentum->Clone("purity_tru0_recF_pBOTHER_momentum");
	TH1F *purity_tru0_recFM_pBOTHER_theta = evt_tru0_recFM_pBOTHER_theta->Clone("purity_tru0_recFM_pBOTHER_theta");
	TH1F *purity_tru0_recFM_pBOTHER_momentum = evt_tru0_recFM_pBOTHER_momentum->Clone("purity_tru0_recFM_pBOTHER_momentum");
	TH1F *purity_tru0_recFMQ_pBOTHER_theta = evt_tru0_recFMQ_pBOTHER_theta->Clone("purity_tru0_recFMQ_pBOTHER_theta");
	TH1F *purity_tru0_recFMQ_pBOTHER_momentum = evt_tru0_recFMQ_pBOTHER_momentum->Clone("purity_tru0_recFMQ_pBOTHER_momentum");
		purity_tru0_recF_pBOTHER_momentum->Divide(evt_tru0_recF_p0_momentum);
		purity_tru0_recF_pBOTHER_theta->Divide(evt_tru0_recF_p0_theta);
		purity_tru0_recFM_pBOTHER_momentum->Divide(evt_tru0_recFM_p0_momentum);
		purity_tru0_recFM_pBOTHER_theta->Divide(evt_tru0_recFM_p0_theta);
		purity_tru0_recFMQ_pBOTHER_momentum->Divide(evt_tru0_recFMQ_p0_momentum);
		purity_tru0_recFMQ_pBOTHER_theta->Divide(evt_tru0_recFMQ_p0_theta);
			purity_tru0_recF_pBOTHER_theta->SetLineWidth(1);
			purity_tru0_recF_pBOTHER_theta->SetLineColor(purity_pBOTHER_color);
			purity_tru0_recF_pBOTHER_theta->SetFillColor(purity_pBOTHER_color);
			purity_tru0_recF_pBOTHER_momentum->SetLineWidth(1);
			purity_tru0_recF_pBOTHER_momentum->SetLineColor(purity_pBOTHER_color);
			purity_tru0_recF_pBOTHER_momentum->SetFillColor(purity_pBOTHER_color);
			purity_tru0_recFM_pBOTHER_theta->SetLineWidth(1);
			purity_tru0_recFM_pBOTHER_theta->SetLineColor(purity_pBOTHER_color);
			purity_tru0_recFM_pBOTHER_theta->SetFillColor(purity_pBOTHER_color);
			purity_tru0_recFM_pBOTHER_momentum->SetLineWidth(1);
			purity_tru0_recFM_pBOTHER_momentum->SetLineColor(purity_pBOTHER_color);
			purity_tru0_recFM_pBOTHER_momentum->SetFillColor(purity_pBOTHER_color);
			purity_tru0_recFMQ_pBOTHER_theta->SetLineWidth(1);
			purity_tru0_recFMQ_pBOTHER_theta->SetLineColor(purity_pBOTHER_color);
			purity_tru0_recFMQ_pBOTHER_theta->SetFillColor(purity_pBOTHER_color);
			purity_tru0_recFMQ_pBOTHER_momentum->SetLineWidth(1);
			purity_tru0_recFMQ_pBOTHER_momentum->SetLineColor(purity_pBOTHER_color);
			purity_tru0_recFMQ_pBOTHER_momentum->SetFillColor(purity_pBOTHER_color);
	}
	// F   purity basic
	{
    THStack purity_basicsummary_F_theta("purity_basicsummary_F_theta","Purities after fiducial volume cut;#theta_{#mu} (degrees);purity");
			purity_basicsummary_F_theta.SetMinimum(0);
			purity_basicsummary_F_theta.Add(purity_tru0_recF_pCC_theta);
			purity_basicsummary_F_theta.Add(purity_tru0_recF_pBKG_theta);
		purity_basicsummary_F_theta.Draw("hist");
 		purity_basicsummary_F_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_basicsummary_F_leg->AddEntry(purity_tru0_recF_pCC_theta,"Charged current","f");
			purity_basicsummary_F_leg->AddEntry(purity_tru0_recF_pBKG_theta,"Background","f");
		purity_basicsummary_F_leg->Draw();
	} gPad->GetCanvas()->Print("purity_basicsummary_F_theta.png");
	{
    THStack purity_basicsummary_F_momentum("purity_basicsummary_F_momentum","Purities after fiducial volume cut;p_{#mu} (GeV/c);purity");
			purity_basicsummary_F_momentum.SetMinimum(0);
			purity_basicsummary_F_momentum.Add(purity_tru0_recF_pCC_momentum);
			purity_basicsummary_F_momentum.Add(purity_tru0_recF_pBKG_momentum);
		purity_basicsummary_F_momentum.Draw("hist");
 		purity_basicsummary_F_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_basicsummary_F_leg->AddEntry(purity_tru0_recF_pCC_momentum,"Charged current","f");
			purity_basicsummary_F_leg->AddEntry(purity_tru0_recF_pBKG_momentum,"Background","f");
		purity_basicsummary_F_leg->Draw();
	} gPad->GetCanvas()->Print("purity_basicsummary_F_momentum.png");
	// FM  purity basic
	{
    THStack purity_basicsummary_FM_theta("purity_basicsummary_FM_theta","Purities after MINOS best match cut;#theta_{#mu} (degrees);purity");
			purity_basicsummary_FM_theta.SetMinimum(0);
			purity_basicsummary_FM_theta.Add(purity_tru0_recFM_pCC_theta);
			purity_basicsummary_FM_theta.Add(purity_tru0_recFM_pBKG_theta);
		purity_basicsummary_FM_theta.Draw("hist");
 		purity_basicsummary_FM_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_basicsummary_FM_leg->AddEntry(purity_tru0_recFM_pCC_theta,"Charged current","f");
			purity_basicsummary_FM_leg->AddEntry(purity_tru0_recFM_pBKG_theta,"Background","f");
		purity_basicsummary_FM_leg->Draw();
	} gPad->GetCanvas()->Print("purity_basicsummary_FM_theta.png");
	{
    THStack purity_basicsummary_FM_momentum("purity_basicsummary_FM_momentum","Purities after MINOS best match cut;p_{#mu} (GeV/c);purity");
			purity_basicsummary_FM_momentum.SetMinimum(0);
			purity_basicsummary_FM_momentum.Add(purity_tru0_recFM_pCC_momentum);
			purity_basicsummary_FM_momentum.Add(purity_tru0_recFM_pBKG_momentum);
		purity_basicsummary_FM_momentum.Draw("hist");
 		purity_basicsummary_FM_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_basicsummary_FM_leg->AddEntry(purity_tru0_recFM_pCC_momentum,"Charged current","f");
			purity_basicsummary_FM_leg->AddEntry(purity_tru0_recFM_pBKG_momentum,"Background","f");
		purity_basicsummary_FM_leg->Draw();
	} gPad->GetCanvas()->Print("purity_basicsummary_FM_momentum.png");
	// FMQ purity basic
	{
    THStack purity_basicsummary_FMQ_theta("purity_basicsummary_FMQ_theta","Purities after all cuts;#theta_{#mu} (degrees);purity");
			purity_basicsummary_FMQ_theta.SetMinimum(0);
			purity_basicsummary_FMQ_theta.Add(purity_tru0_recFMQ_pCC_theta);
			purity_basicsummary_FMQ_theta.Add(purity_tru0_recFMQ_pBKG_theta);
		purity_basicsummary_FMQ_theta.Draw("hist");
 		purity_basicsummary_FMQ_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_basicsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCC_theta,"Charged current","f");
			purity_basicsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pBKG_theta,"Background","f");
		purity_basicsummary_FMQ_leg->Draw();
	} gPad->GetCanvas()->Print("purity_basicsummary_FMQ_theta.png");
	{
    THStack purity_basicsummary_FMQ_momentum("purity_basicsummary_FMQ_momentum","Purities after all cuts;p_{#mu} (GeV/c);purity");
			purity_basicsummary_FMQ_momentum.SetMinimum(0);
			purity_basicsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pCC_momentum);
			purity_basicsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pBKG_momentum);
		purity_basicsummary_FMQ_momentum.Draw("hist");
 		purity_basicsummary_FMQ_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_basicsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCC_momentum,"Charged current","f");
			purity_basicsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pBKG_momentum,"Background","f");
		purity_basicsummary_FMQ_leg->Draw();
	} gPad->GetCanvas()->Print("purity_basicsummary_FMQ_momentum.png");
	// F   purity specific
	{
    THStack purity_detailsummary_F_theta("purity_detailsummary_F_theta","Purities after fiducial volume cut;#theta_{#mu} (degrees);purity");
			purity_detailsummary_F_theta.SetMinimum(0);
			purity_detailsummary_F_theta.Add(purity_tru0_recF_pCCDIS_theta);
			purity_detailsummary_F_theta.Add(purity_tru0_recF_pCCRES_theta);
			purity_detailsummary_F_theta.Add(purity_tru0_recF_pCCQE_theta);
			purity_detailsummary_F_theta.Add(purity_tru0_recF_pCCOTHER_theta);
			purity_detailsummary_F_theta.Add(purity_tru0_recF_pNC_theta);
			purity_detailsummary_F_theta.Add(purity_tru0_recF_pWS_theta);
			purity_detailsummary_F_theta.Add(purity_tru0_recF_pBOTHER_theta);
		purity_detailsummary_F_theta.Draw("hist");
 		purity_detailsummary_F_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pCCDIS_theta,"Deep inelastic","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pCCRES_theta,"Resonant","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pCCQE_theta,"Quasi-elastic","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pCCOTHER_theta,"Other charged","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pNC_theta,"Neutral current","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pWS_theta,"Wrong sign","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pBOTHER_theta,"Other background","f");
    purity_detailsummary_F_leg->Draw();
	} gPad->GetCanvas()->Print("purity_detailsummary_F_theta.png");
	{
    THStack purity_detailsummary_F_momentum("purity_detailsummary_F_momentum","Purities after fiducial volume cut;p_{#mu} (GeV/c);purity");
			purity_detailsummary_F_momentum.SetMinimum(0);
			purity_detailsummary_F_momentum.Add(purity_tru0_recF_pCCDIS_momentum);
			purity_detailsummary_F_momentum.Add(purity_tru0_recF_pCCRES_momentum);
			purity_detailsummary_F_momentum.Add(purity_tru0_recF_pCCQE_momentum);
			purity_detailsummary_F_momentum.Add(purity_tru0_recF_pCCOTHER_momentum);
			purity_detailsummary_F_momentum.Add(purity_tru0_recF_pNC_momentum);
			purity_detailsummary_F_momentum.Add(purity_tru0_recF_pWS_momentum);
			purity_detailsummary_F_momentum.Add(purity_tru0_recF_pBOTHER_momentum);
		purity_detailsummary_F_momentum.Draw("hist");
 		purity_detailsummary_F_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pCCDIS_momentum,"Deep inelastic","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pCCRES_momentum,"Resonant","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pCCQE_momentum,"Quasi-elastic","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pCCOTHER_momentum,"Other charged","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pNC_momentum,"Neutral current","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pWS_momentum,"Wrong sign","f");
			purity_detailsummary_F_leg->AddEntry(purity_tru0_recF_pBOTHER_momentum,"Other background","f");
    purity_detailsummary_F_leg->Draw();
	} gPad->GetCanvas()->Print("purity_detailsummary_F_momentum.png");
	// FM  purity specific
	{
    THStack purity_detailsummary_FM_theta("purity_detailsummary_FM_theta","Purities after MINOS best match cut;#theta_{#mu} (degrees);purity");
			purity_detailsummary_FM_theta.SetMinimum(0);
			purity_detailsummary_FM_theta.Add(purity_tru0_recFM_pCCDIS_theta);
			purity_detailsummary_FM_theta.Add(purity_tru0_recFM_pCCRES_theta);
			purity_detailsummary_FM_theta.Add(purity_tru0_recFM_pCCQE_theta);
			purity_detailsummary_FM_theta.Add(purity_tru0_recFM_pCCOTHER_theta);
			purity_detailsummary_FM_theta.Add(purity_tru0_recFM_pNC_theta);
			purity_detailsummary_FM_theta.Add(purity_tru0_recFM_pWS_theta);
			purity_detailsummary_FM_theta.Add(purity_tru0_recFM_pBOTHER_theta);
		purity_detailsummary_FM_theta.Draw("hist");
 		purity_detailsummary_FM_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pCCDIS_theta,"Deep inelastic","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pCCRES_theta,"Resonant","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pCCQE_theta,"Quasi-elastic","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pCCOTHER_theta,"Other charged","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pNC_theta,"Neutral current","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pWS_theta,"Wrong sign","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pBOTHER_theta,"Other background","f");
    purity_detailsummary_FM_leg->Draw();
	} gPad->GetCanvas()->Print("purity_detailsummary_FM_theta.png");
	{
    THStack purity_detailsummary_FM_momentum("purity_detailsummary_FM_momentum","Purities after MINOS best match cut;p_{#mu} (GeV/c);purity");
			purity_detailsummary_FM_momentum.SetMinimum(0);
			purity_detailsummary_FM_momentum.Add(purity_tru0_recFM_pCCDIS_momentum);
			purity_detailsummary_FM_momentum.Add(purity_tru0_recFM_pCCRES_momentum);
			purity_detailsummary_FM_momentum.Add(purity_tru0_recFM_pCCQE_momentum);
			purity_detailsummary_FM_momentum.Add(purity_tru0_recFM_pCCOTHER_momentum);
			purity_detailsummary_FM_momentum.Add(purity_tru0_recFM_pNC_momentum);
			purity_detailsummary_FM_momentum.Add(purity_tru0_recFM_pWS_momentum);
			purity_detailsummary_FM_momentum.Add(purity_tru0_recFM_pBOTHER_momentum);
		purity_detailsummary_FM_momentum.Draw("hist");
 		purity_detailsummary_FM_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pCCDIS_momentum,"Deep inelastic","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pCCRES_momentum,"Resonant","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pCCQE_momentum,"Quasi-elastic","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pCCOTHER_momentum,"Other charged","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pNC_momentum,"Neutral current","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pWS_momentum,"Wrong sign","f");
			purity_detailsummary_FM_leg->AddEntry(purity_tru0_recFM_pBOTHER_momentum,"Other background","f");
    purity_detailsummary_FM_leg->Draw();
	} gPad->GetCanvas()->Print("purity_detailsummary_FM_momentum.png");
	// FMQ purity specific
	{
    THStack purity_detailsummary_FMQ_theta("purity_detailsummary_FMQ_theta","Purities after all cuts;#theta_{#mu} (degrees);purity");
			purity_detailsummary_FMQ_theta.SetMinimum(0);
			purity_detailsummary_FMQ_theta.Add(purity_tru0_recFMQ_pCCDIS_theta);
			purity_detailsummary_FMQ_theta.Add(purity_tru0_recFMQ_pCCRES_theta);
			purity_detailsummary_FMQ_theta.Add(purity_tru0_recFMQ_pCCQE_theta);
			purity_detailsummary_FMQ_theta.Add(purity_tru0_recFMQ_pCCOTHER_theta);
			purity_detailsummary_FMQ_theta.Add(purity_tru0_recFMQ_pNC_theta);
			purity_detailsummary_FMQ_theta.Add(purity_tru0_recFMQ_pBOTHER_theta);
			purity_detailsummary_FMQ_theta.Add(purity_tru0_recFMQ_pWS_theta);
		purity_detailsummary_FMQ_theta.Draw("hist");
 		purity_detailsummary_FMQ_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCCDIS_theta,"Deep inelastic","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCCRES_theta,"Resonant","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCCQE_theta,"Quasi-elastic","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCCOTHER_theta,"Other charged","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pNC_theta,"Neutral current","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pWS_theta,"Wrong sign","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pBOTHER_theta,"Other background","f");
    purity_detailsummary_FMQ_leg->Draw();
	} gPad->GetCanvas()->Print("purity_detailsummary_FMQ_theta.png");
	{
    THStack purity_detailsummary_FMQ_momentum("purity_detailsummary_FMQ_momentum","Purities after all cuts;p_{#mu} (GeV/c);purity");
			purity_detailsummary_FMQ_momentum.SetMinimum(0);
			purity_detailsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pCCDIS_momentum);
			purity_detailsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pCCRES_momentum);
			purity_detailsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pCCQE_momentum);
			purity_detailsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pCCOTHER_momentum);
			purity_detailsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pNC_momentum);
			purity_detailsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pWS_momentum);
			purity_detailsummary_FMQ_momentum.Add(purity_tru0_recFMQ_pBOTHER_momentum);
		purity_detailsummary_FMQ_momentum.Draw("hist");
 		purity_detailsummary_FMQ_leg = new TLegend(0.6,0.2,0.79,0.39);
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCCDIS_momentum,"Deep inelastic","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCCRES_momentum,"Resonant","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCCQE_momentum,"Quasi-elastic","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pCCOTHER_momentum,"Other charged","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pNC_momentum,"Neutral current","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pWS_momentum,"Wrong sign","f");
			purity_detailsummary_FMQ_leg->AddEntry(purity_tru0_recFMQ_pBOTHER_momentum,"Other background","f");
    purity_detailsummary_FMQ_leg->Draw();
	} gPad->GetCanvas()->Print("purity_detailsummary_FMQ_momentum.png");
	// F   backgrounds
	{
    THStack background_summary_recF_theta("background_summary_recF_theta","Background per total after fiducial volume cut;#theta_{#mu} (degrees);purity");
			background_summary_recF_theta.SetMinimum(0);
			background_summary_recF_theta.Add(purity_tru0_recF_pNC_theta);
			background_summary_recF_theta.Add(purity_tru0_recF_pWS_theta);
			background_summary_recF_theta.Add(purity_tru0_recF_pBOTHER_theta);
		background_summary_recF_theta.Draw("hist");
 		background_summary_recF_leg = new TLegend(0.6,0.2,0.79,0.39);
			background_summary_recF_leg->AddEntry(purity_tru0_recF_pNC_theta,"Neutral current","f");
			background_summary_recF_leg->AddEntry(purity_tru0_recF_pWS_theta,"Wrong sign","f");
			background_summary_recF_leg->AddEntry(purity_tru0_recF_pBOTHER_theta,"Other","f");
		background_summary_recF_leg->Draw();
	} gPad->GetCanvas()->Print("background_summary_recF_theta.png");
	{
    THStack background_summary_recF_momentum("background_summary_recF_momentum","Background per total after fiducial volume cut;p_{#mu} (GeV/c);purity");
			background_summary_recF_momentum.SetMinimum(0);
			background_summary_recF_momentum.Add(purity_tru0_recF_pNC_momentum);
			background_summary_recF_momentum.Add(purity_tru0_recF_pWS_momentum);
			background_summary_recF_momentum.Add(purity_tru0_recF_pBOTHER_momentum);
		background_summary_recF_momentum.Draw("hist");
 		background_summary_recF_leg = new TLegend(0.6,0.2,0.79,0.39);
			background_summary_recF_leg->AddEntry(purity_tru0_recF_pNC_momentum,"Neutral current","f");
			background_summary_recF_leg->AddEntry(purity_tru0_recF_pWS_momentum,"Wrong sign","f");
			background_summary_recF_leg->AddEntry(purity_tru0_recF_pBOTHER_momentum,"Other","f");
		background_summary_recF_leg->Draw();
	} gPad->GetCanvas()->Print("background_summary_recF_momentum.png");
	// FM  backgrounds
	{
    THStack background_summary_recFM_theta("background_summary_recFM_theta","Background per total after MINOS best match cut;#theta_{#mu} (degrees);purity");
			background_summary_recFM_theta.SetMinimum(0);
			background_summary_recFM_theta.Add(purity_tru0_recFM_pNC_theta);
			background_summary_recFM_theta.Add(purity_tru0_recFM_pWS_theta);
			background_summary_recFM_theta.Add(purity_tru0_recFM_pBOTHER_theta);
		background_summary_recFM_theta.Draw("hist");
 		background_summary_recFM_leg = new TLegend(0.6,0.2,0.79,0.39);
			background_summary_recFM_leg->AddEntry(purity_tru0_recFM_pNC_theta,"Neutral current","f");
			background_summary_recFM_leg->AddEntry(purity_tru0_recFM_pWS_theta,"Wrong sign","f");
			background_summary_recFM_leg->AddEntry(purity_tru0_recFM_pBOTHER_theta,"Other","f");
		background_summary_recFM_leg->Draw();
	} gPad->GetCanvas()->Print("background_summary_recFM_theta.png");
	{
    THStack background_summary_recFM_momentum("background_summary_recFM_momentum","Background per total after MINOS best match cut;p_{#mu} (GeV/c);purity");
			background_summary_recFM_momentum.SetMinimum(0);
			background_summary_recFM_momentum.Add(purity_tru0_recFM_pNC_momentum);
			background_summary_recFM_momentum.Add(purity_tru0_recFM_pWS_momentum);
			background_summary_recFM_momentum.Add(purity_tru0_recFM_pBOTHER_momentum);
		background_summary_recFM_momentum.Draw("hist");
 		background_summary_recFM_leg = new TLegend(0.6,0.2,0.79,0.39);
			background_summary_recFM_leg->AddEntry(purity_tru0_recFM_pNC_momentum,"Neutral current","f");
			background_summary_recFM_leg->AddEntry(purity_tru0_recFM_pWS_momentum,"Wrong sign","f");
			background_summary_recFM_leg->AddEntry(purity_tru0_recFM_pBOTHER_momentum,"Other","f");
		background_summary_recFM_leg->Draw();
	} gPad->GetCanvas()->Print("background_summary_recFM_momentum.png");
	// FMQ backgrounds
	{
    THStack background_summary_recFMQ_theta("background_summary_recFMQ_theta","Background per total after all cuts;#theta_{#mu} (degrees);purity");
			background_summary_recFMQ_theta.SetMinimum(0);
			background_summary_recFMQ_theta.Add(purity_tru0_recFMQ_pNC_theta);
			background_summary_recFMQ_theta.Add(purity_tru0_recFMQ_pWS_theta);
			background_summary_recFMQ_theta.Add(purity_tru0_recFMQ_pBOTHER_theta);
		background_summary_recFMQ_theta.Draw("hist");
 		background_summary_recFMQ_leg = new TLegend(0.6,0.2,0.79,0.39);
			background_summary_recFMQ_leg->AddEntry(purity_tru0_recFMQ_pNC_theta,"Neutral current","f");
			background_summary_recFMQ_leg->AddEntry(purity_tru0_recFMQ_pWS_theta,"Wrong sign","f");
			background_summary_recFMQ_leg->AddEntry(purity_tru0_recFMQ_pBOTHER_theta,"Other","f");
		background_summary_recFMQ_leg->Draw();
	} gPad->GetCanvas()->Print("background_summary_recFMQ_theta.png");
	{
    THStack background_summary_recFMQ_momentum("background_summary_recFMQ_momentum","Background per total after all cuts;p_{#mu} (GeV/c);purity");
			background_summary_recFMQ_momentum.SetMinimum(0);
			background_summary_recFMQ_momentum.Add(purity_tru0_recFMQ_pNC_momentum);
			background_summary_recFMQ_momentum.Add(purity_tru0_recFMQ_pWS_momentum);
			background_summary_recFMQ_momentum.Add(purity_tru0_recFMQ_pBOTHER_momentum);
		background_summary_recFMQ_momentum.Draw("hist");
 		background_summary_recFMQ_leg = new TLegend(0.6,0.2,0.79,0.39);
			background_summary_recFMQ_leg->AddEntry(purity_tru0_recFMQ_pNC_momentum,"Neutral current","f");
			background_summary_recFMQ_leg->AddEntry(purity_tru0_recFMQ_pWS_momentum,"Wrong sign","f");
			background_summary_recFMQ_leg->AddEntry(purity_tru0_recFMQ_pBOTHER_momentum,"Other","f");
		background_summary_recFMQ_leg->Draw();
	} gPad->GetCanvas()->Print("background_summary_recFMQ_momentum.png");

	// counts
	{
	int count_p0_color       = 39;
	int count_pCC_color      = 31;
	int count_pCCQE_color    = 30;
	int count_pCCRES_color   = 31;
	int count_pCCDIS_color   = 33;
	int count_pCCOTHER_color = 32;
	int count_pNC_color      = 20;
	int count_pWS_color      = 24;
	int count_pBOTHER_color  = 27;

	evt_tru0_recF_p0_theta->SetLineColor(count_p0_color);
	evt_tru0_recF_p0_momentum->SetLineColor(count_p0_color);
	evt_tru0_recFM_p0_theta->SetLineColor(count_p0_color);
	evt_tru0_recFM_p0_momentum->SetLineColor(count_p0_color);
	evt_tru0_recFMQ_p0_theta->SetLineColor(count_p0_color);
	evt_tru0_recFMQ_p0_momentum->SetLineColor(count_p0_color);

	evt_tru0_recF_pCC_theta->SetLineColor(count_pCC_color);
	evt_tru0_recF_pCC_momentum->SetLineColor(count_pCC_color);
	evt_tru0_recFM_pCC_theta->SetLineColor(count_pCC_color);
	evt_tru0_recFM_pCC_momentum->SetLineColor(count_pCC_color);
	evt_tru0_recFMQ_pCC_theta->SetLineColor(count_pCC_color);
	evt_tru0_recFMQ_pCC_momentum->SetLineColor(count_pCC_color);

	evt_tru0_recF_pCCRES_theta->SetLineColor(count_pCCRES_color);
	evt_tru0_recF_pCCRES_momentum->SetLineColor(count_pCCRES_color);
	evt_tru0_recFM_pCCRES_theta->SetLineColor(count_pCCRES_color);
	evt_tru0_recFM_pCCRES_momentum->SetLineColor(count_pCCRES_color);
	evt_tru0_recFMQ_pCCRES_theta->SetLineColor(count_pCCRES_color);
	evt_tru0_recFMQ_pCCRES_momentum->SetLineColor(count_pCCRES_color);

	evt_tru0_recF_pCCDIS_theta->SetLineColor(count_pCCDIS_color);
	evt_tru0_recF_pCCDIS_momentum->SetLineColor(count_pCCDIS_color);
	evt_tru0_recFM_pCCDIS_theta->SetLineColor(count_pCCDIS_color);
	evt_tru0_recFM_pCCDIS_momentum->SetLineColor(count_pCCDIS_color);
	evt_tru0_recFMQ_pCCDIS_theta->SetLineColor(count_pCCDIS_color);
	evt_tru0_recFMQ_pCCDIS_momentum->SetLineColor(count_pCCDIS_color);

	evt_tru0_recF_pCCQE_theta->SetLineColor(count_pCCQE_color);
	evt_tru0_recF_pCCQE_momentum->SetLineColor(count_pCCQE_color);
	evt_tru0_recFM_pCCQE_theta->SetLineColor(count_pCCQE_color);
	evt_tru0_recFM_pCCQE_momentum->SetLineColor(count_pCCQE_color);
	evt_tru0_recFMQ_pCCQE_theta->SetLineColor(count_pCCQE_color);
	evt_tru0_recFMQ_pCCQE_momentum->SetLineColor(count_pCCQE_color);
	
	evt_tru0_recF_pCCOTHER_theta->SetLineColor(count_pCCOTHER_color);
	evt_tru0_recF_pCCOTHER_momentum->SetLineColor(count_pCCOTHER_color);
	evt_tru0_recFM_pCCOTHER_theta->SetLineColor(count_pCCOTHER_color);
	evt_tru0_recFM_pCCOTHER_momentum->SetLineColor(count_pCCOTHER_color);
	evt_tru0_recFMQ_pCCOTHER_theta->SetLineColor(count_pCCOTHER_color);
	evt_tru0_recFMQ_pCCOTHER_momentum->SetLineColor(count_pCCOTHER_color);

	evt_tru0_recF_pNC_theta->SetLineColor(count_pNC_color);
	evt_tru0_recF_pNC_momentum->SetLineColor(count_pNC_color);
	evt_tru0_recFM_pNC_theta->SetLineColor(count_pNC_color);
	evt_tru0_recFM_pNC_momentum->SetLineColor(count_pNC_color);
	evt_tru0_recFMQ_pNC_theta->SetLineColor(count_pNC_color);
	evt_tru0_recFMQ_pNC_momentum->SetLineColor(count_pNC_color);

	evt_tru0_recF_pWS_theta->SetLineColor(count_pWS_color);
	evt_tru0_recF_pWS_momentum->SetLineColor(count_pWS_color);
	evt_tru0_recFM_pWS_theta->SetLineColor(count_pWS_color);
	evt_tru0_recFM_pWS_momentum->SetLineColor(count_pWS_color);
	evt_tru0_recFMQ_pWS_theta->SetLineColor(count_pWS_color);
	evt_tru0_recFMQ_pWS_momentum->SetLineColor(count_pWS_color);

	evt_tru0_recF_pBOTHER_theta->SetLineColor(count_pBOTHER_color);
	evt_tru0_recF_pBOTHER_momentum->SetLineColor(count_pBOTHER_color);
	evt_tru0_recFM_pBOTHER_theta->SetLineColor(count_pBOTHER_color);
	evt_tru0_recFM_pBOTHER_momentum->SetLineColor(count_pBOTHER_color);
	evt_tru0_recFMQ_pBOTHER_theta->SetLineColor(count_pBOTHER_color);
	evt_tru0_recFMQ_pBOTHER_momentum->SetLineColor(count_pBOTHER_color);
	}
	// F   counts
	{
    THStack count_recF_theta("count_recF_theta","Events after fiducial volume cut;#theta_{#mu} (degrees)");
		count_recF_theta.Add(evt_tru0_recF_p0_theta);
		count_recF_theta.Add(evt_tru0_recF_pCC_theta);
 		count_recF_theta.Add(evt_tru0_recF_pNC_theta);
		count_recF_theta.Add(evt_tru0_recF_pWS_theta);
		count_recF_theta.Add(evt_tru0_recF_pBOTHER_theta);
		count_recF_theta.Draw("hist nostack");
		leg = new TLegend(0.6,0.6,0.8,0.8);
			leg->AddEntry(evt_tru0_recF_p0_theta,"All","l");
			leg->AddEntry(evt_tru0_recF_pCC_theta,"Charged current","l");
			leg->AddEntry(evt_tru0_recF_pNC_theta,"Neutral current","l");
			leg->AddEntry(evt_tru0_recF_pWS_theta,"Wrong sign","l");
			leg->AddEntry(evt_tru0_recF_pBOTHER_theta,"Other","l");
    leg->Draw();
	} gPad->GetCanvas()->Print("count_recF_theta.png");
	{
    THStack count_recF_momentum("count_recF_momentum","Events after fiducial volume cut;p_{#mu} (GeV/c)");
		count_recF_momentum.Add(evt_tru0_recF_p0_momentum);
		count_recF_momentum.Add(evt_tru0_recF_pCC_momentum);
 		count_recF_momentum.Add(evt_tru0_recF_pNC_momentum);
		count_recF_momentum.Add(evt_tru0_recF_pWS_momentum);
		count_recF_momentum.Add(evt_tru0_recF_pBOTHER_momentum);
		count_recF_momentum.Draw("hist nostack");
		leg = new TLegend(0.6,0.6,0.8,0.8);
			leg->AddEntry(evt_tru0_recF_p0_momentum,"All","l");
			leg->AddEntry(evt_tru0_recF_pCC_momentum,"Charged current","l");
			leg->AddEntry(evt_tru0_recF_pNC_momentum,"Neutral current","l");
			leg->AddEntry(evt_tru0_recF_pWS_momentum,"Wrong sign","l");
			leg->AddEntry(evt_tru0_recF_pBOTHER_momentum,"Other","l");
    leg->Draw();
	} gPad->GetCanvas()->Print("count_recF_momentum.png");
	// FM  counts
	{
    THStack count_recFM_theta("count_recFM_theta","Events after MINOS best match cut;#theta_{#mu} (degrees)");
		count_recFM_theta.Add(evt_tru0_recFM_p0_theta);
		count_recFM_theta.Add(evt_tru0_recFM_pCC_theta);
 		count_recFM_theta.Add(evt_tru0_recFM_pNC_theta);
		count_recFM_theta.Add(evt_tru0_recFM_pWS_theta);
		count_recFM_theta.Add(evt_tru0_recFM_pBOTHER_theta);
		count_recFM_theta.Draw("hist nostack");
		leg = new TLegend(0.6,0.6,0.8,0.8);
			leg->AddEntry(evt_tru0_recFM_p0_theta,"All","l");
			leg->AddEntry(evt_tru0_recFM_pCC_theta,"Charged current","l");
			leg->AddEntry(evt_tru0_recFM_pNC_theta,"Neutral current","l");
			leg->AddEntry(evt_tru0_recFM_pWS_theta,"Wrong sign","l");
			leg->AddEntry(evt_tru0_recFM_pBOTHER_theta,"Other","l");
    leg->Draw();
	} gPad->GetCanvas()->Print("count_recFM_theta.png");
	{
    THStack count_recFM_momentum("count_recFM_momentum","Events after MINOS best match cut;p_{#mu} (GeV/c)");
		count_recFM_momentum.Add(evt_tru0_recFM_p0_momentum);
		count_recFM_momentum.Add(evt_tru0_recFM_pCC_momentum);
 		count_recFM_momentum.Add(evt_tru0_recFM_pNC_momentum);
		count_recFM_momentum.Add(evt_tru0_recFM_pWS_momentum);
		count_recFM_momentum.Add(evt_tru0_recFM_pBOTHER_momentum);
		count_recFM_momentum.Draw("hist nostack");
		leg = new TLegend(0.6,0.6,0.8,0.8);
			leg->AddEntry(evt_tru0_recFM_p0_momentum,"All","l");
			leg->AddEntry(evt_tru0_recFM_pCC_momentum,"Charged current","l");
			leg->AddEntry(evt_tru0_recFM_pNC_momentum,"Neutral current","l");
			leg->AddEntry(evt_tru0_recFM_pWS_momentum,"Wrong sign","l");
			leg->AddEntry(evt_tru0_recFM_pBOTHER_momentum,"Other","l");
    leg->Draw();
	} gPad->GetCanvas()->Print("count_recFM_momentum.png");
	// FMQ counts
	{
    THStack count_recFMQ_theta("count_recFMQ_theta","Events after MINOS charge cut;#theta_{#mu} (degrees)");
		count_recFMQ_theta.Add(evt_tru0_recFMQ_p0_theta);
		count_recFMQ_theta.Add(evt_tru0_recFMQ_pCC_theta);
 		count_recFMQ_theta.Add(evt_tru0_recFMQ_pNC_theta);
		count_recFMQ_theta.Add(evt_tru0_recFMQ_pWS_theta);
		count_recFMQ_theta.Add(evt_tru0_recFMQ_pBOTHER_theta);
		count_recFMQ_theta.Draw("hist nostack");
		leg = new TLegend(0.6,0.6,0.8,0.8);
			leg->AddEntry(evt_tru0_recFMQ_p0_theta,"All","l");
			leg->AddEntry(evt_tru0_recFMQ_pCC_theta,"Charged current","l");
			leg->AddEntry(evt_tru0_recFMQ_pNC_theta,"Neutral current","l");
			leg->AddEntry(evt_tru0_recFMQ_pWS_theta,"Wrong sign","l");
			leg->AddEntry(evt_tru0_recFMQ_pBOTHER_theta,"Other","l");
    leg->Draw();
	} gPad->GetCanvas()->Print("count_recFMQ_theta.png");
	{
    THStack count_recFMQ_momentum("count_recFMQ_momentum","Events after MINOS charge cut;p_{#mu} (GeV/c)");
		count_recFMQ_momentum.Add(evt_tru0_recFMQ_p0_momentum);
		count_recFMQ_momentum.Add(evt_tru0_recFMQ_pCC_momentum);
 		count_recFMQ_momentum.Add(evt_tru0_recFMQ_pNC_momentum);
		count_recFMQ_momentum.Add(evt_tru0_recFMQ_pWS_momentum);
		count_recFMQ_momentum.Add(evt_tru0_recFMQ_pBOTHER_momentum);
		count_recFMQ_momentum.Draw("hist nostack");
		leg = new TLegend(0.6,0.6,0.8,0.8);
			leg->AddEntry(evt_tru0_recFMQ_p0_momentum,"All","l");
			leg->AddEntry(evt_tru0_recFMQ_pCC_momentum,"Charged current","l");
			leg->AddEntry(evt_tru0_recFMQ_pNC_momentum,"Neutral current","l");
			leg->AddEntry(evt_tru0_recFMQ_pWS_momentum,"Wrong sign","l");
			leg->AddEntry(evt_tru0_recFMQ_pBOTHER_momentum,"Other","l");
    leg->Draw();
	} gPad->GetCanvas()->Print("count_recFMQ_momentum.png");

  /*
	// all F counts
	new TCanvas;	{
		evt_tru0_recF_p0_theta->SetTitle("Events after fiducial volume cuts");
		evt_tru0_recF_p0_theta->SetXTitle("#theta_{#mu, true} (GeV/c)");
		evt_tru0_recF_p0_theta->SetYTitle("events"); 
		evt_tru0_recF_p0_theta->SetLineColor(33);
		evt_tru0_recF_p0_theta->Draw("hist");
		evt_tru0_recF_pCC_theta->SetLineColor(42);
		evt_tru0_recF_pCC_theta->Draw("hist same");
		evt_tru0_recF_pCCDIS_theta->SetLineColor(42);
		evt_tru0_recF_pCCDIS_theta->Draw("hist same");
		evt_tru0_recF_pCCRES_theta->SetLineColor(30);
		evt_tru0_recF_pCCRES_theta->Draw("hist same");
		evt_tru0_recF_pCCDIS_theta->SetLineColor(46);
		evt_tru0_recF_pCCDIS_theta->Draw("hist same");
		evt_tru0_recF_pCCQE_theta->SetLineColor(38);
		evt_tru0_recF_pCCQE_theta->Draw("hist same");
		evt_tru0_recF_pNC_theta->SetLineColor(27);
		evt_tru0_recF_pNC_theta->Draw("hist same");
		evt_tru0_recF_pBOTHER_theta->SetLineColor(47);
		evt_tru0_recF_pBOTHER_theta->Draw("hist same");
	} gPad->GetCanvas()->Print("count_recoF_theta.png");
	new TCanvas;	{
		evt_tru0_recF_p0_momentum->SetTitle("Events after fiducial volume cuts");
		evt_tru0_recF_p0_momentum->SetXTitle("p_{#mu, true} (GeV/c)");
		evt_tru0_recF_p0_momentum->SetYTitle("events");
		evt_tru0_recF_p0_momentum->SetLineColor(33);
		evt_tru0_recF_p0_momentum->Draw("hist");
		evt_tru0_recF_pCC_momentum->SetLineColor(42);
		evt_tru0_recF_pCC_momentum->Draw("hist same");
		evt_tru0_recF_pCCDIS_momentum->SetLineColor(42);
		evt_tru0_recF_pCCDIS_momentum->Draw("hist same");
		evt_tru0_recF_pCCRES_momentum->SetLineColor(30);
		evt_tru0_recF_pCCRES_momentum->Draw("hist same");
		evt_tru0_recF_pCCDIS_momentum->SetLineColor(46);
		evt_tru0_recF_pCCDIS_momentum->Draw("hist same");
		evt_tru0_recF_pCCQE_momentum->SetLineColor(38);
		evt_tru0_recF_pCCQE_momentum->Draw("hist same");
		evt_tru0_recF_pNC_momentum->SetLineColor(27);
		evt_tru0_recF_pNC_momentum->Draw("hist same");
		evt_tru0_recF_pBOTHER_momentum->SetLineColor(47);
		evt_tru0_recF_pBOTHER_momentum->Draw("hist same");
	} gPad->GetCanvas()->Print("count_recoF_momentum.png");
	// all FM counts
	new TCanvas;	{
		gStyle->SetOptStat(0);
		evt_tru0_recFM_p0_theta->SetTitle("Events after matching cuts");
		evt_tru0_recFM_p0_theta->SetXTitle("#theta_{#mu, true} (GeV/c)");
		evt_tru0_recFM_p0_theta->SetYTitle("events"); 
		evt_tru0_recFM_p0_theta->SetLineColor(33);
		evt_tru0_recFM_p0_theta->Draw("hist");
		evt_tru0_recFM_pCC_theta->SetLineColor(42);
		evt_tru0_recFM_pCC_theta->Draw("hist same");
		evt_tru0_recFM_pCCDIS_theta->SetLineColor(42);
		evt_tru0_recFM_pCCDIS_theta->Draw("hist same");
		evt_tru0_recFM_pCCRES_theta->SetLineColor(30);
		evt_tru0_recFM_pCCRES_theta->Draw("hist same");
		evt_tru0_recFM_pCCDIS_theta->SetLineColor(46);
		evt_tru0_recFM_pCCDIS_theta->Draw("hist same");
		evt_tru0_recFM_pCCQE_theta->SetLineColor(38);
		evt_tru0_recFM_pCCQE_theta->Draw("hist same");
		evt_tru0_recFM_pNC_theta->SetLineColor(27);
		evt_tru0_recFM_pNC_theta->Draw("hist same");
		evt_tru0_recFM_pBOTHER_theta->SetLineColor(47);
		evt_tru0_recFM_pBOTHER_theta->Draw("hist same");
	} gPad->GetCanvas()->Print("count_recoFM_theta.png");
	new TCanvas;	{
		gStyle->SetOptStat(0);
		evt_tru0_recFM_p0_momentum->SetTitle("Events after matching cuts");
		evt_tru0_recFM_p0_momentum->SetXTitle("p_{#mu, true} (GeV/c)");
		evt_tru0_recFM_p0_momentum->SetYTitle("events");
		evt_tru0_recFM_p0_momentum->SetLineColor(33);
		evt_tru0_recFM_p0_momentum->Draw("hist");
		evt_tru0_recFM_pCC_momentum->SetLineColor(42);
		evt_tru0_recFM_pCC_momentum->Draw("hist same");
		evt_tru0_recFM_pCCDIS_momentum->SetLineColor(42);
		evt_tru0_recFM_pCCDIS_momentum->Draw("hist same");
		evt_tru0_recFM_pCCRES_momentum->SetLineColor(30);
		evt_tru0_recFM_pCCRES_momentum->Draw("hist same");
		evt_tru0_recFM_pCCDIS_momentum->SetLineColor(46);
		evt_tru0_recFM_pCCDIS_momentum->Draw("hist same");
		evt_tru0_recFM_pCCQE_momentum->SetLineColor(38);
		evt_tru0_recFM_pCCQE_momentum->Draw("hist same");
		evt_tru0_recFM_pNC_momentum->SetLineColor(27);
		evt_tru0_recFM_pNC_momentum->Draw("hist same");
		evt_tru0_recFM_pBOTHER_momentum->SetLineColor(47);
		evt_tru0_recFM_pBOTHER_momentum->Draw("hist same");
	} gPad->GetCanvas()->Print("count_recoFM_momentum.png");
	// all FMQ counts (final)
	new TCanvas;	{
		gStyle->SetOptStat(0);
		evt_tru0_recFMQ_p0_theta->SetTitle("Events after all cuts");
		evt_tru0_recFMQ_p0_theta->SetXTitle("#theta_{#mu, true} (GeV/c)");
		evt_tru0_recFMQ_p0_theta->SetYTitle("events"); 
		evt_tru0_recFMQ_p0_theta->SetLineColor(33);
		evt_tru0_recFMQ_p0_theta->Draw("hist");
		evt_tru0_recFMQ_pCC_theta->SetLineColor(42);
		evt_tru0_recFMQ_pCC_theta->Draw("hist same");
		evt_tru0_recFMQ_pCCDIS_theta->SetLineColor(42);
		evt_tru0_recFMQ_pCCDIS_theta->Draw("hist same");
		evt_tru0_recFMQ_pCCRES_theta->SetLineColor(30);
		evt_tru0_recFMQ_pCCRES_theta->Draw("hist same");
		evt_tru0_recFMQ_pCCDIS_theta->SetLineColor(46);
		evt_tru0_recFMQ_pCCDIS_theta->Draw("hist same");
		evt_tru0_recFMQ_pCCQE_theta->SetLineColor(38);
		evt_tru0_recFMQ_pCCQE_theta->Draw("hist same");
		evt_tru0_recFMQ_pNC_theta->SetLineColor(27);
		evt_tru0_recFMQ_pNC_theta->Draw("hist same");
		evt_tru0_recFMQ_pBOTHER_theta->SetLineColor(47);
		evt_tru0_recFMQ_pBOTHER_theta->Draw("hist same");
	} gPad->GetCanvas()->Print("count_recoFMQ_theta.png");
	new TCanvas;	{
		gStyle->SetOptStat(0);
		evt_tru0_recFMQ_p0_momentum->SetTitle("Events after all cuts");
		evt_tru0_recFMQ_p0_momentum->SetXTitle("p_{#mu, true} (GeV/c)");
		evt_tru0_recFMQ_p0_momentum->SetYTitle("events");
		evt_tru0_recFMQ_p0_momentum->SetLineColor(33);
		evt_tru0_recFMQ_p0_momentum->Draw("hist");
		evt_tru0_recFMQ_pCC_momentum->SetLineColor(42);
		evt_tru0_recFMQ_pCC_momentum->Draw("hist same");
		evt_tru0_recFMQ_pCCDIS_momentum->SetLineColor(42);
		evt_tru0_recFMQ_pCCDIS_momentum->Draw("hist same");
		evt_tru0_recFMQ_pCCRES_momentum->SetLineColor(30);
		evt_tru0_recFMQ_pCCRES_momentum->Draw("hist same");
		evt_tru0_recFMQ_pCCDIS_momentum->SetLineColor(46);
		evt_tru0_recFMQ_pCCDIS_momentum->Draw("hist same");
		evt_tru0_recFMQ_pCCQE_momentum->SetLineColor(38);
		evt_tru0_recFMQ_pCCQE_momentum->Draw("hist same");
		evt_tru0_recFMQ_pNC_momentum->SetLineColor(27);
		evt_tru0_recFMQ_pNC_momentum->Draw("hist same");
		evt_tru0_recFMQ_pBOTHER_momentum->SetLineColor(47);
		evt_tru0_recFMQ_pBOTHER_momentum->Draw("hist same");
	} gPad->GetCanvas()->Print("count_recoFMQ_momentum.png");
  */
  /*
	// all kinds, bkgd purities
	new TCanvas;	{
		purity_tru0_recFMQ_pCC_theta->SetTitle("Purities after all cuts");
		purity_tru0_recFMQ_pCC_theta->SetXTitle("#theta_{#mu, true} (GeV/c)");
		purity_tru0_recFMQ_pCC_theta->SetYTitle("events");
		purity_tru0_recFMQ_pCC_theta->SetMaximum(1);
		purity_tru0_recFMQ_pCC_theta->SetMinimum(0);
		purity_tru0_recFMQ_pCC_theta->SetStats(0);
		purity_tru0_recFMQ_pCC_theta->SetLineColor(33);
		purity_tru0_recFMQ_pCC_theta->Draw("hist");
		purity_tru0_recFMQ_pCC_theta->SetLineColor(42);
		purity_tru0_recFMQ_pCC_theta->Draw("hist same");
		purity_tru0_recFMQ_pCCDIS_theta->SetLineColor(42);
		purity_tru0_recFMQ_pCCDIS_theta->Draw("hist same");
		purity_tru0_recFMQ_pCCRES_theta->SetLineColor(30);
		purity_tru0_recFMQ_pCCRES_theta->Draw("hist same");
		purity_tru0_recFMQ_pCCDIS_theta->SetLineColor(46);
		purity_tru0_recFMQ_pCCDIS_theta->Draw("hist same");
		purity_tru0_recFMQ_pCCQE_theta->SetLineColor(38);
		purity_tru0_recFMQ_pCCQE_theta->Draw("hist same");
		purity_tru0_recFMQ_pNC_theta->SetLineColor(27);
		purity_tru0_recFMQ_pNC_theta->Draw("hist same");
		purity_tru0_recFMQ_pBOTHER_theta->SetLineColor(47);
		purity_tru0_recFMQ_pBOTHER_theta->Draw("hist same");
	} gPad->GetCanvas()->Print("purities_all_theta.png");
	new TCanvas;	{
		gStyle->SetOptStat(0);
		purity_tru0_recFMQ_pCC_momentum->SetTitle("Purities after all cuts");
		purity_tru0_recFMQ_pCC_momentum->SetXTitle("p_{#mu, true} (GeV/c)");
		purity_tru0_recFMQ_pCC_momentum->SetYTitle("events");
		purity_tru0_recFMQ_pCC_momentum->SetMaximum(1);
		purity_tru0_recFMQ_pCC_momentum->SetMinimum(0);
		purity_tru0_recFMQ_pCC_momentum->SetStats(0);
		purity_tru0_recFMQ_pCC_momentum->SetLineColor(33);
		purity_tru0_recFMQ_pCC_momentum->Draw("hist");
		purity_tru0_recFMQ_pCC_momentum->SetLineColor(42);
		purity_tru0_recFMQ_pCC_momentum->Draw("hist same");
		purity_tru0_recFMQ_pCCDIS_momentum->SetLineColor(42);
		purity_tru0_recFMQ_pCCDIS_momentum->Draw("hist same");
		purity_tru0_recFMQ_pCCRES_momentum->SetLineColor(30);
		purity_tru0_recFMQ_pCCRES_momentum->Draw("hist same");
		purity_tru0_recFMQ_pCCDIS_momentum->SetLineColor(46);
		purity_tru0_recFMQ_pCCDIS_momentum->Draw("hist same");
		purity_tru0_recFMQ_pCCQE_momentum->SetLineColor(38);
		purity_tru0_recFMQ_pCCQE_momentum->Draw("hist same");
		purity_tru0_recFMQ_pNC_momentum->SetLineColor(27);
		purity_tru0_recFMQ_pNC_momentum->Draw("hist same");
		purity_tru0_recFMQ_pBOTHER_momentum->SetLineColor(47);
		purity_tru0_recFMQ_pBOTHER_momentum->Draw("hist same");
	} gPad->GetCanvas()->Print("purities_all_momentum.png");
  */
  /*
	// reco FMQ counts
	new TCanvas;	{
		evt_tru0_recFMQ_p0_theta->Divide(evt_truF_rec0_pCC_theta); // use same denominator as ArgoNeuT alone
		evt_tru0_recFMQ_p0_theta->SetTitle("[NEW] CC #nu_{#mu} muon ArgoNeuT+MINOS match and q reconstruction probability (after cuts)");
		evt_tru0_recFMQ_p0_theta->SetXTitle("#theta_{#mu, true} (degrees)");
		evt_tru0_recFMQ_p0_theta->SetYTitle("probability");
		evt_tru0_recFMQ_p0_theta->SetMaximum(1);
		evt_tru0_recFMQ_p0_theta->SetMinimum(0);
		evt_tru0_recFMQ_p0_theta->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_p0_theta.png");
	new TCanvas;	{
		evt_tru0_recFMQ_p0_momentum->Divide(evt_truF_rec0_pCC_momentum); // use same denominator as ArgoNeuT alone
		evt_tru0_recFMQ_p0_momentum->SetTitle("[NEW] CC #nu_{#mu} muon ArgoNeuT+MINOS match and q reconstruction probability (after cuts)");
		evt_tru0_recFMQ_p0_momentum->SetXTitle("P_{#mu, true} (degrees)");
		evt_tru0_recFMQ_p0_momentum->SetYTitle("probability");
		evt_tru0_recFMQ_p0_momentum->SetMaximum(1);
		evt_tru0_recFMQ_p0_momentum->SetMinimum(0);
		evt_tru0_recFMQ_p0_momentum->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_p0_momentum.png");
  // CC FMQ counts
	new TCanvas;	{
		evt_tru0_recFMQ_pCC_theta->SetTitle("Charged current events after all cuts");
		evt_tru0_recFMQ_pCC_theta->SetXTitle("#theta_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pCC_theta->SetYTitle("events");
		evt_tru0_recFMQ_pCC_theta->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pCC_theta.png");
	new TCanvas;	{
		evt_tru0_recFMQ_pCC_momentum->SetTitle("Charged current events after all cuts");
		evt_tru0_recFMQ_pCC_momentum->SetXTitle("P_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pCC_momentum->SetYTitle("events");
		evt_tru0_recFMQ_pCC_momentum->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pCC_momentum.png");
  // CCQE FMQ counts
	new TCanvas;	{
		evt_tru0_recFMQ_pCCQE_theta->SetTitle("Charged current events after all cuts");
		evt_tru0_recFMQ_pCCQE_theta->SetXTitle("#theta_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pCCQE_theta->SetYTitle("events");
		evt_tru0_recFMQ_pCCQE_theta->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pCCQE_theta.png");
	new TCanvas;	{
		evt_tru0_recFMQ_pCCQE_momentum->SetTitle("Charged current QE events after all cuts");
		evt_tru0_recFMQ_pCCQE_momentum->SetXTitle("P_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pCCQE_momentum->SetYTitle("events");
		evt_tru0_recFMQ_pCCQE_momentum->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pCCQE_momentum.png");
  // CCRES FMQ counts
	new TCanvas;	{
		evt_tru0_recFMQ_pCCRES_theta->SetTitle("Charged current resonant events after all cuts");
		evt_tru0_recFMQ_pCCRES_theta->SetXTitle("#theta_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pCCRES_theta->SetYTitle("events");
		evt_tru0_recFMQ_pCCRES_theta->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pCCRES_theta.png");
	new TCanvas;	{
		evt_tru0_recFMQ_pCCRES_momentum->SetTitle("Charged current resonant events after all cuts");
		evt_tru0_recFMQ_pCCRES_momentum->SetXTitle("P_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pCCRES_momentum->SetYTitle("events");
		evt_tru0_recFMQ_pCCRES_momentum->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pCCRES_momentum.png");
	// CCDIS FMQ counts
	new TCanvas;	{
		evt_tru0_recFMQ_pCCDIS_theta->SetTitle("Charged current DIS events after all cuts");
		evt_tru0_recFMQ_pCCDIS_theta->SetXTitle("#theta_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pCCDIS_theta->SetYTitle("events");
		evt_tru0_recFMQ_pCCDIS_theta->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pCCDIS_theta.png");
	new TCanvas;	{
		evt_tru0_recFMQ_pCCDIS_momentum->SetTitle("Charged current DIS events after all cuts");
		evt_tru0_recFMQ_pCCDIS_momentum->SetXTitle("P_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pCCDIS_momentum->SetYTitle("events");
		evt_tru0_recFMQ_pCCDIS_momentum->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pCCDIS_momentum.png");
	// NC FMQ counts
	new TCanvas;	{
		evt_tru0_recFMQ_pNC_theta->SetTitle("Neutral current events after all cuts");
		evt_tru0_recFMQ_pNC_theta->SetXTitle("#theta_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pNC_theta->SetYTitle("events");
		evt_tru0_recFMQ_pNC_theta->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pNC_theta.png");
	new TCanvas;	{
		evt_tru0_recFMQ_pNC_momentum->SetTitle("Neutral current events after all cuts");
		evt_tru0_recFMQ_pNC_momentum->SetXTitle("P_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pNC_momentum->SetYTitle("events");
		evt_tru0_recFMQ_pNC_momentum->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pNC_momentum.png");
	// OTHER FMQ counts
	new TCanvas;	{
		evt_tru0_recFMQ_pBOTHER_theta->SetTitle("Other events after all cuts");
		evt_tru0_recFMQ_pBOTHER_theta->SetXTitle("#theta_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pBOTHER_theta->SetYTitle("events");
		evt_tru0_recFMQ_pBOTHER_theta->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pBOTHER_theta.png");
	new TCanvas;	{
		evt_tru0_recFMQ_pBOTHER_momentum->SetTitle("Other events after all cuts");
		evt_tru0_recFMQ_pBOTHER_momentum->SetXTitle("P_{#mu, true} (degrees)");
		evt_tru0_recFMQ_pBOTHER_momentum->SetYTitle("events");
		evt_tru0_recFMQ_pBOTHER_momentum->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFMQ_pBOTHER_momentum.png");

	// efficiency of neutrino in fiducial volume being matched in MINOS
	new TCanvas;	{
		evt_tru0_recFM_p0_theta->Divide(evt_truF_rec0_pCC_theta); // use same denominator as ArgoNeuT alone
		evt_tru0_recFM_p0_theta->SetTitle("[NEW] CC #nu_{#mu} muon ArgoNeuT+MINOS match reconstruction probability (after cuts)");
		evt_tru0_recFM_p0_theta->SetXTitle("#theta_{#mu, true} (degrees)");
		evt_tru0_recFM_p0_theta->SetYTitle("probability");
		evt_tru0_recFM_p0_theta->SetMaximum(1);
		evt_tru0_recFM_p0_theta->SetMinimum(0);
		evt_tru0_recFM_p0_theta->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFM_p0_theta.png");
	new TCanvas;	{
		evt_tru0_recFM_p0_momentum->Divide(evt_truF_rec0_pCC_momentum); // use same denominator as ArgoNeuT alone
		evt_tru0_recFM_p0_momentum->SetTitle("[NEW] CC #nu_{#mu} muon ArgoNeuT+MINOS match reconstruction probability (after cuts)");
		evt_tru0_recFM_p0_momentum->SetXTitle("P_{#mu, true} (degrees)");
		evt_tru0_recFM_p0_momentum->SetYTitle("probability");
		evt_tru0_recFM_p0_momentum->SetMaximum(1);
		evt_tru0_recFM_p0_momentum->SetMinimum(0);
		evt_tru0_recFM_p0_momentum->DrawCopy("E1");
	} gPad->GetCanvas()->Print("evt_tru0_recFM_p0_momentum.png");

	new TCanvas;	{
		evt_truFQ_rec0_p0_theta->Divide(evt_truF_rec0_pCC_theta);
		evt_truFQ_rec0_p0_theta->SetTitle("CC #nu_{#mu} muon MINOS reconstruction probability (after cuts)");
		evt_truFQ_rec0_p0_theta->SetXTitle("#theta_{#mu, true}");
		evt_truFQ_rec0_p0_theta->SetYTitle("probability");
		evt_truFQ_rec0_p0_theta->SetMaximum(1);
		evt_truFQ_rec0_p0_theta->SetMinimum(0);
		evt_truFQ_rec0_p0_theta->Draw();
	} gPad->GetCanvas()->Print("evt_truFQ_rec0_p0_theta.png");
	new TCanvas;	{
		evt_truFQ_rec0_p0_momentum->Divide(evt_truF_rec0_pCC_momentum);
		evt_truFQ_rec0_p0_momentum->SetTitle("CC #nu_{#mu} muon MINOS reconstruction probability (after cuts)");
		evt_truFQ_rec0_p0_momentum->SetXTitle("P_{#mu} (GeV/c)");
		evt_truFQ_rec0_p0_momentum->SetYTitle("probability");
		evt_truFQ_rec0_p0_momentum->SetMaximum(1);
		evt_truFQ_rec0_p0_momentum->SetMinimum(0);
		evt_truFQ_rec0_p0_momentum->Draw();
	} gPad->GetCanvas()->Print("evt_truFQ_rec0_p0_momentum.png"); 



	// ratio of events reconstructed in fiducial volume to to true events in fiducial volume 
	new TCanvas;	{
		TH1F *pArgoFV_reco_theta_probability = evt_tru0_recF_p0_theta->Clone("pArgoFV_reco_theta_probability");
		pArgoFV_reco_theta_probability->Divide(evt_truF_rec0_pCC_theta); // use same denominator as ArgoNeuT alone
		pArgoFV_reco_theta_probability->SetTitle("Ratio of reconstructed events in fiducial volume to true events in fiducial volume");
		pArgoFV_reco_theta_probability->SetXTitle("#theta_{#mu, true} (degrees)");
		pArgoFV_reco_theta_probability->SetYTitle("ratio");
		pArgoFV_reco_theta_probability->SetMaximum(2);
		pArgoFV_reco_theta_probability->SetMinimum(0);
		pArgoFV_reco_theta_probability->DrawCopy("E1");
	} gPad->GetCanvas()->Print("ratio_recF-truF_theta.png");
	new TCanvas;	{
		TH1F *pArgoFV_reco_momentum_probability = evt_tru0_recF_p0_momentum->Clone("pArgoFV_reco_momentum_probability");
		pArgoFV_reco_momentum_probability->Divide(evt_truF_rec0_pCC_momentum); // use same denominator as ArgoNeuT alone
		pArgoFV_reco_momentum_probability->SetTitle("Ratio of reconstructed events in fiducial volume to true events in fiducial volume");
		pArgoFV_reco_momentum_probability->SetXTitle("P_{#mu, true} (degrees)");
		pArgoFV_reco_momentum_probability->SetYTitle("radio");
		pArgoFV_reco_momentum_probability->SetMaximum(2);
		pArgoFV_reco_momentum_probability->SetMinimum(0);
		pArgoFV_reco_momentum_probability->DrawCopy("E1");
	} gPad->GetCanvas()->Print("ratio_recF-truF_momentum.png");

	// ratios of events reconstructed in fiducial volume to all true events reconstructed in fiducial volume
	new TCanvas;	{
		gStyle->SetOptStat(0);
		TH1F *pArgoFV_reco_theta_efficiency = evt_tru0_recF_p0_theta->Clone("pArgoFV_reco_theta_efficiency");
		TH1F *pArgoFV_MINOSm_reco_theta_efficiency = evt_tru0_recFM_p0_theta->Clone("pArgoFV_MINOSm_reco_theta_efficiency");
		TH1F *pArgoFV_MINOSmq_reco_theta_efficiency = evt_tru0_recFMQ_p0_theta->Clone("pArgoFV_MINOSmq_reco_theta_efficiency");
		pArgoFV_MINOSm_reco_theta_efficiency->Divide(evt_truF_recF_p0_theta); // use same denominator as ArgoNeuT alone
		pArgoFV_MINOSmq_reco_theta_efficiency->Divide(evt_truF_recF_p0_theta); // use same denominator as ArgoNeuT alone
		pArgoFV_reco_theta_efficiency->Divide(evt_truF_recF_p0_theta); // use same denominator as ArgoNeuT alone
		pArgoFV_reco_theta_efficiency->SetTitle("reco / recotrue for each cut"); 	// not working...
		pArgoFV_reco_theta_efficiency->SetXTitle("#theta_{#mu, true} (GeV/c)"); 									// not working...
		pArgoFV_reco_theta_efficiency->SetYTitle("ratio"); // not working...
		pArgoFV_reco_theta_efficiency->SetMinimum(0);
		pArgoFV_reco_theta_efficiency->SetFillColor(33);
		pArgoFV_reco_theta_efficiency->SetLineColor(33);
		pArgoFV_reco_theta_efficiency->Draw();
		pArgoFV_MINOSm_reco_theta_efficiency->SetFillColor(42);
		pArgoFV_MINOSm_reco_theta_efficiency->SetLineColor(42);
		pArgoFV_MINOSm_reco_theta_efficiency->Draw("same");
		pArgoFV_MINOSmq_reco_theta_efficiency->SetFillColor(30);
		pArgoFV_MINOSmq_reco_theta_efficiency->SetLineColor(30);
		pArgoFV_MINOSmq_reco_theta_efficiency->Draw("same");
		leg->Draw();
	} gPad->GetCanvas()->Print("ratio_recX-rectruF_theta.png");
	new TCanvas;	{
		gStyle->SetOptStat(0);
		TH1F *pArgoFV_reco_momentum_efficiency = evt_tru0_recF_p0_momentum->Clone("pArgoFV_reco_momentum_efficiency");
		TH1F *pArgoFV_MINOSm_reco_momentum_efficiency = evt_tru0_recFM_p0_momentum->Clone("pArgoFV_MINOSm_reco_momentum_efficiency");
		TH1F *pArgoFV_MINOSmq_reco_momentum_efficiency = evt_tru0_recFMQ_p0_momentum->Clone("pArgoFV_MINOSmq_reco_momentum_efficiency");
		pArgoFV_MINOSm_reco_momentum_efficiency->Divide(evt_truF_recF_p0_momentum); // use same denominator as ArgoNeuT alone
		pArgoFV_MINOSmq_reco_momentum_efficiency->Divide(evt_truF_recF_p0_momentum); // use same denominator as ArgoNeuT alone
		pArgoFV_reco_momentum_efficiency->Divide(evt_truF_recF_p0_momentum); // use same denominator as ArgoNeuT alone
		pArgoFV_reco_momentum_efficiency->SetTitle("reco / recotrue for each cut"); 	// not working...
		pArgoFV_reco_momentum_efficiency->SetXTitle("P_{#mu, true} (GeV/c)"); 												// not working...
		pArgoFV_reco_momentum_efficiency->SetYTitle("ratio"); 	// not working...
		pArgoFV_reco_momentum_efficiency->SetMinimum(0);
		pArgoFV_reco_momentum_efficiency->SetFillColor(33);
		pArgoFV_reco_momentum_efficiency->SetLineColor(33);
		pArgoFV_reco_momentum_efficiency->Draw();
		pArgoFV_MINOSm_reco_momentum_efficiency->SetFillColor(42);
		pArgoFV_MINOSm_reco_momentum_efficiency->SetLineColor(42);
		pArgoFV_MINOSm_reco_momentum_efficiency->Draw("same");
		pArgoFV_MINOSmq_reco_momentum_efficiency->SetFillColor(30);
		pArgoFV_MINOSmq_reco_momentum_efficiency->SetLineColor(30);
		pArgoFV_MINOSmq_reco_momentum_efficiency->SetStats(0);
		pArgoFV_MINOSmq_reco_momentum_efficiency->Draw("same");
		leg->Draw();
	} gPad->GetCanvas()->Print("ratio_recX-rectruF_momentum.png");



	*/

}

// project ArgoNeuT track into MINOS and calculate difference
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

// set style of plots
void set_plot_style()
{
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	
	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);
}

// method to give weight of event
double Weigh( double enu_true , double weight[14]) { 
	if( enu_true >=  0 && enu_true <  3 )
		return weight[0];
	if( enu_true >=  3 && enu_true <  4 )
		return weight[1];
	if( enu_true >=  4 && enu_true <  5 )
		return weight[2];
	if( enu_true >=  5 && enu_true <  7 )
		return weight[3];
	if( enu_true >=  7 && enu_true <  9 )
		return weight[4];
	if( enu_true >=  9 && enu_true < 12 )
		return weight[5];
	if( enu_true >= 12 && enu_true < 15 )
		return weight[6];
	if( enu_true >= 15 && enu_true < 18 )
		return weight[7];
	if( enu_true >= 18 && enu_true < 22 )
		return weight[8];
	if( enu_true >= 22 && enu_true < 26 )
		return weight[9];
	if( enu_true >= 26 && enu_true < 30 )
		return weight[10];
	if( enu_true >= 30 && enu_true < 36 )
		return weight[11];
	if( enu_true >= 36 && enu_true < 42 )
		return weight[12];
	if( enu_true >= 42 && enu_true < 50 )
		return weight[13];
}
