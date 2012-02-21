#define datakin_cxx
#include "datakin.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "cutvariables.h"
void datakin::Loop()
{
 gROOT->SetStyle("Plain");
   gStyle->SetPalette(1,0); 
   gStyle->SetLineWidth(2); 
   gStyle->SetHistLineWidth(3);
   gStyle->SetOptFit(kTRUE);
  // gStyle->SetOptStat(kFALSE);
//    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptStat(1000001110);
   gStyle->SetOptFit(0011);
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(.032,"Z");
    gStyle->SetStatFontSize(.02);
    gStyle->SetTitleFontSize(.044);
   TH1::AddDirectory(false);
   set_plot_style();
   
   TH2D* study = new TH2D("", "",100,0.5,1.,100,0,50);
   TH2D* study_qe = new TH2D("", "",100,0.5,1.,100,0,50);
   TH2D* study_res = new TH2D("", "",100,0.5,1.,100,0,50);
   TH2D* study_dis = new TH2D("", "",100,0.5,1.,100,0,50);
    TH2D* histunfold = new TH2D("nominal", "nominal",20,0.,40,20,0.,40);
    TH2D* hist3unfoldup = new TH2D("up", "up",20,0.,40,20,0.,40);
    TH2D* hist3unfolddown = new TH2D("down", "down",20,0.,40,20,0.,40);
    
    TH1D* histunfold2 = new TH1D("", "",200,-10,10);
    TH2D* histunfold3 = new TH2D("", "",100,60.,120.,100,60.,120.);
    TH1D* histunfold4 = new TH1D("", "",100,-60,60);
    TH2D* histunfold5 = new TH2D("", "",100,60.,120.,100,60.,120.);
    TH1D* histunfold6 = new TH1D("", "",100,-60,60);
    TH2D* histunfold33 = new TH2D("", "",100,60.,120.,100,60.,120.);
    TH1D* histunfold44 = new TH1D("", "",100,-60,60);
    TH2D* histunfold55 = new TH2D("", "",100,60.,120.,100,60.,120.);
    TH1D* histunfold66 = new TH1D("", "",100,-60,60);
//     TH1D* histtest = new TH1D("", "",60,-1,1);
//     TH1D* histtest2 = new TH1D("", "",60,-1,1);
//     TH1D* histtest3 = new TH1D("", "",60,-3.14,3.14);
//     TH2D* histtest4 = new TH2D("", "",100,-1,1,100,-1,1);
    
    TH2D* histunfoldmom = new TH2D("nominal", "nominal",25,0.,31.25,25,0.,31.25);
    TH2D* hist3unfoldmueup = new TH2D("up", "up",25,0.,31.25,25,0.,31.25);
    TH2D* hist3unfoldmuedown = new TH2D("down", "down",25,0.,31.25,25,0.,31.25);
    
    TH1D* histunfoldmom2 = new TH1D("", "",100,-5,5);
    TH1D* histunfoldmom22 = new TH1D("", "",100,-5,5);
    TH2D* histmomsyst = new TH2D("", "",20,0.,25,50,0,2);
    TH2D* histmomsyst2 = new TH2D("", "",20,0.,25,50,-.5,.5);
    TH2D* histmomsyst3 = new TH2D("", "",20,0.,25,50,-.2,.2);

    TH2D* histthetasyst = new TH2D("", "",18,0.,36,50,-5.,5.);
    TH2D* histthetasyst2 = new TH2D("", "",18,0.,36,50,-1.,1.);
    TH2D* histthetasyst3 = new TH2D("", "",18,0.,36,50,-1.,1.);
    
    TH2D* histmomsyst2ext = new TH2D("", "",40,0.,50,50,-.5,.5);
    TH2D* histthetasyst2ext = new TH2D("", "",36,0.,72,50,-1.,1.);
    
    TH1F* vtxhist = new TH1F("", "",80,-10,10);
    TH1F* vtxhist2 = new TH1F("", "",80,-10,10);
    TH1F* vtxhist3 = new TH1F("", "",80,-10,10);
    TH1F* vtxhist4 = new TH1F("", "",160,-20,20);
    
    TH2F* vtxhistt = new TH2F("", "",47,0,47,47,0,47);
    TH2F* vtxhistt2 = new TH2F("", "",40,-20,20,40,-20,20);
    TH2F* vtxhistt3 = new TH2F("", "",90,0,90,90,0,90);

    
    TH1D* histeff1 = new TH1D("", "",18,0.,36);
    histeff1->Sumw2();
    TH1D* histeff2 = new TH1D("", "",18,0.,36);
    histeff2->Sumw2();
    TH1D* histeff11 = new TH1D("", "",20,0.,25);
    histeff11->Sumw2();
    TH1D* histeff22 = new TH1D("", "",20,0.,25);
    histeff22->Sumw2();
    
    TH1D* histeff222 = new TH1D("", "",18,0.,36);
    histeff222->Sumw2();
    
    TH1D* histeff2222 = new TH1D("", "",20,0.,25);
    histeff2222->Sumw2();
    
    TH1D* histeff111 = new TH1D("", "",18,0.,36);
    histeff111->Sumw2();
    TH1D* histeff1111 = new TH1D("", "",20,0.,25);
     histeff1111->Sumw2();
    TH1D* histpur1 = new TH1D("", "",18,0.,36);
    histpur1->Sumw2();
    TH1D* histpur2 = new TH1D("", "",18,0.,36);
    histpur2->Sumw2();
    
    TH1D* hit3 = new TH1D("", "",20,0.,25);
    hit3->Sumw2();
    TH1D* hit3iso = new TH1D("isomom", "",20,0.,25);
    hit3iso->Sumw2();
    TH1D* hit4 = new TH1D("", "",18,0.,36);
    hit4->Sumw2();
    TH1D* hit4iso = new TH1D("isocos", "",18,0.,36);
     hit4iso->Sumw2();
    TH1D*  histeff_qe1 = new TH1D("", "",20,0.,25);
    TH1D*  histeff_qe2 = new TH1D("", "",20,0.,25);
    TH1D*  histeff_res1 = new TH1D("", "",20,0.,25);
    TH1D*  histeff_res2 = new TH1D("", "",20,0.,25);
    TH1D*  histeff_dis1 = new TH1D("", "",20,0.,25);
    TH1D*  histeff_dis2 = new TH1D("", "",20,0.,25);
    
    TH1D* histbkgd = new TH1D("", "",18,0.,36);
     histbkgd->Sumw2();
    TH1D* histbkgd2 = new TH1D("", "",20,0.,25);
     histbkgd2->Sumw2();
    
//    

TH1D* hist1 = new TH1D("", "",60,0,30);
TH1D* hist2 = new TH1D("", "",60,0,30);
TH1D* hist3 = new TH1D("", "",60,0,30);
TH1D* hist4 = new TH1D("", "",60,0,30);
TH1D* hist5 = new TH1D("", "",60,0,30);
TH1D* hist6 = new TH1D("", "",60,0,30);

TH1D* hit2 = new TH1D("", "",50,0,50);
TH1D* histunweighted = new TH1D("", "",200,0,50);
 hit2->Sumw2();
TH1D* hit2iso = new TH1D("", "",50,0,50);
hit2iso->Sumw2();
TH1D* hist11 = new TH1D("", "",60,0,30);
TH1D* hist22 = new TH1D("", "",60,0,30);
TH1D* hist33 = new TH1D("", "",60,0,30);
TH1D* hist44 = new TH1D("", "",60,0,30);
TH1D* hist55 = new TH1D("", "",60,0,30);
TH1D* hist66 = new TH1D("", "",60,0,30);

TH1D* hist1nt= new TH1D("", "",10,0,10);
TH1D* hist2nt = new TH1D("", "",10,0,10);
TH1D* hist3nt = new TH1D("", "",10,0,10);
TH1D* hist4nt = new TH1D("", "",10,0,10);
TH1D* hist5nt = new TH1D("", "",10,0,10);
TH1D* hist6nt = new TH1D("", "",10,0,10);

TH1D* hist1nvt= new TH1D("", "",10,0,10);
TH1D* hist2nvt = new TH1D("", "",10,0,10);
TH1D* hist3nvt = new TH1D("", "",10,0,10);
TH1D* hist4nvt = new TH1D("", "",10,0,10);
TH1D* hist5nvt = new TH1D("", "",10,0,10);
TH1D* hist6nvt = new TH1D("", "",10,0,10);

TH1D* hist1ntexit= new TH1D("", "",10,0,10);
TH1D* hist2ntexit = new TH1D("", "",10,0,10);
TH1D* hist3ntexit = new TH1D("", "",10,0,10);
TH1D* hist4ntexit = new TH1D("", "",10,0,10);
TH1D* hist5ntexit = new TH1D("", "",10,0,10);
TH1D* hist6ntexit = new TH1D("", "",10,0,10);

TH1D* histdistance = new TH1D("", "",200,0,200);

TH1D* hist1nvtxclusu= new TH1D("", "",15,0,15);
TH1D* hist2nvtxclusu = new TH1D("", "",15,0,15);
TH1D* hist3nvtxclusu = new TH1D("", "",15,0,15);
TH1D* hist4nvtxclusu = new TH1D("", "",15,0,15);
TH1D* hist5nvtxclusu = new TH1D("", "",15,0,15);
TH1D* hist6nvtxclusu = new TH1D("", "",15,0,15);

TH1D* hist1nvtxclusv= new TH1D("", "",15,0,15);
TH1D* hist2nvtxclusv = new TH1D("", "",15,0,15);
TH1D* hist3nvtxclusv = new TH1D("", "",15,0,15);
TH1D* hist4nvtxclusv = new TH1D("", "",15,0,15);
TH1D* hist5nvtxclusv = new TH1D("", "",15,0,15);
TH1D* hist6nvtxclusv = new TH1D("", "",15,0,15);
 
TH1D* hist1nclusu= new TH1D("", "",15,0,15);
TH1D* hist2nclusu = new TH1D("", "",15,0,15);
TH1D* hist3nclusu = new TH1D("", "",15,0,15);
TH1D* hist4nclusu = new TH1D("", "",15,0,15);
TH1D* hist5nclusu = new TH1D("", "",15,0,15);
TH1D* hist6nclusu = new TH1D("", "",15,0,15);

TH1D* hist1nclusv= new TH1D("", "",15,0,15);
TH1D* hist2nclusv = new TH1D("", "",15,0,15);
TH1D* hist3nclusv = new TH1D("", "",15,0,15);
TH1D* hist4nclusv = new TH1D("", "",15,0,15);
TH1D* hist5nclusv = new TH1D("", "",15,0,15);
TH1D* hist6nclusv = new TH1D("", "",15,0,15);
 
TChain ch("KinTree");


//default
// ch.Add("simkinreco_hist_golden1_1.root/analysistree/anatree");
//   ch.Add("simkinreco_hist_golden2_1.root/analysistree/anatree");
//   ch.Add("simkinreco_hist_golden3_1.root/analysistree/anatree");
//   ch.Add("simkinreco_hist_golden4_1.root/analysistree/anatree");
  
  
//w/ mitch spacepoints 
ch.Add("simkinreco_hist_golden1_2.root/analysistree/anatree");  
ch.Add("simkinreco_hist_golden2_2.root/analysistree/anatree"); 
ch.Add("simkinreco_hist_golden3_2.root/analysistree/anatree"); 
ch.Add("simkinreco_hist_golden4_2.root/analysistree/anatree"); 
 
//  ch.Add("simkinreco_hist_iso.root/analysistree/anatree");
//  ch.Add("simkinreco_hist_iso_2.root/analysistree/anatree");
//  ch.Add("simkinreco_hist_iso_3.root/analysistree/anatree");
//  ch.Add("simkinreco_hist_iso_4.root/analysistree/anatree"); 
//   
//    ch.Add("simkinreco_hist_test.root/analysistree/anatree");
  
 datakin aEvent(KinTree);

   Long64_t nentries = KinTree.GetEntries();
   TFile *outfile = new TFile("out.root","RECREATE");
  cout << "Number of Entries " << nentries << endl;
  float theta=0.;
  double numevents=0;
  double numtotalevents=0;
  double numncevents=0;
  double numccevents_beforecuts=0;
  double numccevents_aftercuts=0;
  
  double numccevents_beforecuts_mom=0;
  double numccevents_aftercuts_mom=0;
   double numccevents_beforecuts_theta=0;
  double numccevents_aftercuts_theta=0;
  
  double numccqeevents_beforecuts=0;
  double numccqeevents_aftercuts=0;
  
  double numccresevents_beforecuts=0;
  double numccresevents_aftercuts=0;
  
  double numccdisevents_beforecuts=0;
  double numccdisevents_aftercuts=0;
  double numccevents_beforecuts_recminos=0;
  double numccqeevents_beforecuts_recminos=0;
  double numccresevents_beforecuts_recminos=0;
  double numccdisevents_beforecuts_recminos=0;
  
  double numccevents_aftercuts_recminos=0;
  //obvious cuts
  double numevents_passmatchandcharge=0;
  //not obvious cuts
  double numevents_passntracks=0;
  double numevents_passnclusters=0;
    double totpot;

   TFile *fi=new TFile("profsyst.root","READ");
   TH1F *profthetasyst2 = (TH1F*)fi->Get("histthetasystfinal");
   TH1F *profmomsyst2 = (TH1F*)fi->Get("histmomsystfinal");
   

   double thetasyst[30]= {profthetasyst2->GetBinContent(1), profthetasyst2->GetBinContent(2), profthetasyst2->GetBinContent(3), profthetasyst2->GetBinContent(4), profthetasyst2->GetBinContent(5), profthetasyst2->GetBinContent(6), profthetasyst2->GetBinContent(7), profthetasyst2->GetBinContent(8), profthetasyst2->GetBinContent(9), profthetasyst2->GetBinContent(10), profthetasyst2->GetBinContent(11), profthetasyst2->GetBinContent(12), profthetasyst2->GetBinContent(13), profthetasyst2->GetBinContent(14),profthetasyst2->GetBinContent(15),profthetasyst2->GetBinContent(16),profthetasyst2->GetBinContent(17),profthetasyst2->GetBinContent(18),profthetasyst2->GetBinContent(19),profthetasyst2->GetBinContent(20),profthetasyst2->GetBinContent(21),profthetasyst2->GetBinContent(22),profthetasyst2->GetBinContent(23),profthetasyst2->GetBinContent(24),profthetasyst2->GetBinContent(25),profthetasyst2->GetBinContent(26),profthetasyst2->GetBinContent(27),profthetasyst2->GetBinContent(28),profthetasyst2->GetBinContent(29),profthetasyst2->GetBinContent(30)};
  
   double momsyst[30]= {profmomsyst2->GetBinContent(1), profmomsyst2->GetBinContent(2), profmomsyst2->GetBinContent(3), profmomsyst2->GetBinContent(4), profmomsyst2->GetBinContent(5), profmomsyst2->GetBinContent(6), profmomsyst2->GetBinContent(7), profmomsyst2->GetBinContent(8), profmomsyst2->GetBinContent(9), profmomsyst2->GetBinContent(10), profmomsyst2->GetBinContent(11), profmomsyst2->GetBinContent(12), profmomsyst2->GetBinContent(13), profmomsyst2->GetBinContent(14),profmomsyst2->GetBinContent(15),profmomsyst2->GetBinContent(16),profmomsyst2->GetBinContent(17),profmomsyst2->GetBinContent(18),profmomsyst2->GetBinContent(19),profmomsyst2->GetBinContent(20),profmomsyst2->GetBinContent(21),profmomsyst2->GetBinContent(22),profmomsyst2->GetBinContent(23),profmomsyst2->GetBinContent(24),profmomsyst2->GetBinContent(25),profmomsyst2->GetBinContent(26),profmomsyst2->GetBinContent(27),profmomsyst2->GetBinContent(28),profmomsyst2->GetBinContent(29),profmomsyst2->GetBinContent(30)};

    
    
   TFile *f=new TFile("energylost.root","READ");
   TH1F *energylost = (TH1F*)f->Get("energylost");

   double energyloss[14]= {energylost->GetBinContent(1), energylost->GetBinContent(2), energylost->GetBinContent(3), energylost->GetBinContent(4), energylost->GetBinContent(5), energylost->GetBinContent(6), energylost->GetBinContent(7), energylost->GetBinContent(8), energylost->GetBinContent(9), energylost->GetBinContent(10), energylost->GetBinContent(11), energylost->GetBinContent(12), energylost->GetBinContent(13), energylost->GetBinContent(14)}; 
    
    TFile *f=new TFile("numu_numode_final.root","READ");
   TH1F *reweight = (TH1F*)f->Get("histdiv");

   double weight[14]= {reweight->GetBinContent(1), reweight->GetBinContent(2), reweight->GetBinContent(3), reweight->GetBinContent(4), reweight->GetBinContent(5), reweight->GetBinContent(6), reweight->GetBinContent(7), reweight->GetBinContent(8), reweight->GetBinContent(9), reweight->GetBinContent(10), reweight->GetBinContent(11), reweight->GetBinContent(12), reweight->GetBinContent(13), reweight->GetBinContent(14)};
   
double rweight=0.;  
    
for (int event = 0; event<nentries; event++) 
{
if(event%10000 == 0)
cout<<event<<" / "<<nentries<<endl;
aEvent.GetEntry(event);  
totpot+=aEvent.pot/100.;


double dcosy=(aEvent.trackstart_dcosy_reco*TMath::Cos(.0583497))+(aEvent.trackstart_dcosz_reco*TMath::Sin(.0583497));
double dcosz=(aEvent.trackstart_dcosz_reco*TMath::Cos(.0583497))-(aEvent.trackstart_dcosy_reco*TMath::Sin(.0583497));

double dcosytruth=(aEvent.lep_dcosy_truth*TMath::Cos(.0583497))+(aEvent.lep_dcosz_truth*TMath::Sin(.0583497));
double dcosztruth=(aEvent.lep_dcosz_truth*TMath::Cos(.0583497))-(aEvent.lep_dcosy_truth*TMath::Sin(.0583497));


if(aEvent.enu_truth>=0&&aEvent.enu_truth<3)
rweight=weight[0];
if(aEvent.enu_truth>=3&&aEvent.enu_truth<4)
rweight=weight[1];
if(aEvent.enu_truth>=4&&aEvent.enu_truth<5)
rweight=weight[2];
if(aEvent.enu_truth>=5&&aEvent.enu_truth<7)
rweight=weight[3];
if(aEvent.enu_truth>=7&&aEvent.enu_truth<9)
rweight=weight[4];
if(aEvent.enu_truth>=9&&aEvent.enu_truth<12)
rweight=weight[5];
if(aEvent.enu_truth>=12&&aEvent.enu_truth<15)
rweight=weight[6];
if(aEvent.enu_truth>=15&&aEvent.enu_truth<18)
rweight=weight[7];
if(aEvent.enu_truth>=18&&aEvent.enu_truth<22)
rweight=weight[8];
if(aEvent.enu_truth>=22&&aEvent.enu_truth<26)
rweight=weight[9];
if(aEvent.enu_truth>=26&&aEvent.enu_truth<30)
rweight=weight[10];
if(aEvent.enu_truth>=30&&aEvent.enu_truth<36)
rweight=weight[11];
if(aEvent.enu_truth>=36&&aEvent.enu_truth<42)
rweight=weight[12];
if(aEvent.enu_truth>=42&&aEvent.enu_truth<50)
rweight=weight[13];

double isoweight=1.;
if(aEvent.hitnuc_truth==2212)
isoweight=22./20.;
if(aEvent.hitnuc_truth==2112)
isoweight=18./20.;




if(
aEvent.nuvtxx_truth>0.+FIDX_CUT 
&& aEvent.nuvtxx_truth<47.-FIDX_CUT 
&& aEvent.nuvtxy_truth>-20.+FIDY_CUT 
&& aEvent.nuvtxy_truth<20.-FIDY_CUT 
&& aEvent.nuvtxz_truth>0.+FIDZup_CUT 
&& aEvent.nuvtxz_truth<90.-FIDZdown_CUT 
)
{

numtotalevents=numtotalevents+(1.*rweight);

hist1->Fill(aEvent.enu_truth,rweight);
if(aEvent.ccnc_truth==1)
{
hist2->Fill(aEvent.enu_truth,rweight);
numncevents=numncevents+(1.*rweight);
}

if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
{
histeff1->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
histeff11->Fill(aEvent.lep_mom_truth,rweight);
hist3->Fill(aEvent.enu_truth,rweight);
hit2->Fill(aEvent.enu_truth,rweight);
histunweighted->Fill(aEvent.enu_truth);
hit2iso->Fill(aEvent.enu_truth,rweight*isoweight);

hit3->Fill(aEvent.lep_mom_truth,rweight);
hit3iso->Fill(aEvent.lep_mom_truth,rweight*isoweight);

hit4->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
hit4iso->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight*isoweight);


numccevents_beforecuts=numccevents_beforecuts+(1.*rweight);


if(aEvent.lep_mom_truth<25.)
numccevents_beforecuts_mom=numccevents_beforecuts_mom+(1.*rweight);

if((180./3.14159)*TMath::ACos(dcosztruth)<36.)
numccevents_beforecuts_theta=numccevents_beforecuts_theta+(1.*rweight);




study->Fill((180./3.14159)*TMath::ACos(dcosztruth),aEvent.lep_mom_truth,rweight);

 if(aEvent.test_charge_minos==-1)
 {
 numccevents_beforecuts_recminos=numccevents_beforecuts_recminos+(1.*rweight);
 histeff111->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
 histeff1111->Fill(aEvent.lep_mom_truth,rweight);
 }
}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
{
hist4->Fill(aEvent.enu_truth,rweight);
numccqeevents_beforecuts=numccqeevents_beforecuts+(1.*rweight);
study_qe->Fill((180./3.14159)*TMath::ACos(dcosztruth),aEvent.lep_mom_truth,rweight);
histeff_qe1->Fill(aEvent.lep_mom_truth,rweight);

if(aEvent.test_charge_minos==-1)
numccqeevents_beforecuts_recminos=numccqeevents_beforecuts_recminos+(1.*rweight);

}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
{
study_res->Fill((180./3.14159)*TMath::ACos(dcosztruth),aEvent.lep_mom_truth,rweight);
histeff_res1->Fill(aEvent.lep_mom_truth,rweight);

hist5->Fill(aEvent.enu_truth,rweight);
numccresevents_beforecuts=numccresevents_beforecuts+(1.*rweight);
if(aEvent.test_charge_minos==-1)
numccresevents_beforecuts_recminos=numccresevents_beforecuts_recminos+(1.*rweight);
}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
{
study_dis->Fill((180./3.14159)*TMath::ACos(dcosztruth),aEvent.lep_mom_truth,rweight);
histeff_dis1->Fill(aEvent.lep_mom_truth,rweight);

hist6->Fill(aEvent.enu_truth,rweight);
numccdisevents_beforecuts=numccdisevents_beforecuts+(1.*rweight);
if(aEvent.test_charge_minos==-1)
numccdisevents_beforecuts_recminos=numccdisevents_beforecuts_recminos+(1.*rweight);
}
}


if(
aEvent.vtxx_reco>0.+FIDX_CUT && aEvent.trackstart_x_reco>0.+FIDX_CUT 
&& aEvent.vtxx_reco<47.-FIDX_CUT && aEvent.trackstart_x_reco<47.-FIDX_CUT  
&& aEvent.vtxy_reco>-20.+FIDY_CUT && aEvent.trackstart_y_reco>-20.+FIDY_CUT 
&& aEvent.vtxy_reco<20.-FIDY_CUT  && aEvent.trackstart_y_reco<20.-FIDY_CUT
&& aEvent.vtxz_reco>0.+FIDZup_CUT && aEvent.trackstart_z_reco>0.+FIDZup_CUT 
&& aEvent.vtxz_reco<90.-FIDZdown_CUT && aEvent.trackstart_z_reco<90.-FIDZdown_CUT 
)
{
double diff_dcosx=aEvent.trackexit_dcosx_reco-aEvent.trk_dcosx_minos;
double diff_dcosy=aEvent.trackexit_dcosy_reco-aEvent.trk_dcosy_minos;
double diff_dcosz=aEvent.trackexit_dcosz_reco-aEvent.trk_dcosz_minos;

//obvious cuts
if(aEvent.nmatched_reco==1&&aEvent.trk_charge_minos<0&&aEvent.trk_mom_minos>0)
{

double lardirectionEnd[3]={aEvent.trackexit_dcosx_reco, aEvent.trackexit_dcosy_reco, aEvent.trackexit_dcosz_reco};
double larEnd[3]= { aEvent.trackexit_x_reco, aEvent.trackexit_y_reco, aEvent.trackexit_z_reco};
double minosdirectionStart[3]={ aEvent.trk_dcosx_minos, aEvent.trk_dcosy_minos, aEvent.trk_dcosz_minos };

double minosStart[3]={aEvent.trk_vtxx_minos, aEvent.trk_vtxy_minos, aEvent.trk_vtxz_minos };
double xdiff,ydiff,rdiff,totaldiff,thetadiff;


project(lardirectionEnd, larEnd, minosdirectionStart,minosStart,xdiff,ydiff,rdiff,totaldiff,thetadiff);

double minosvtx[3] = {(aEvent.trk_vtxx_minos*100.)-117.4, (aEvent.trk_vtxy_minos*100.)+19.3, (aEvent.trk_vtxz_minos*100.)+147.1};

double distance=sqrt(pow(aEvent.vtxx_reco-minosvtx[0],2)+pow(aEvent.vtxy_reco-minosvtx[1],2)+pow(aEvent.vtxz_reco-minosvtx[2],2));




if(rdiff>RDIFF_CUT||thetadiff>THETADIFF_CUT)
continue;

// if(aEvent.mc_pdg_minos!=-99999)
// numccevents_aftercuts_recminos++;
histdistance->Fill(distance,rweight);
MOM_OFFSET=(0.002198*distance)-.04404; 
//       x is distance (cm), y is energy lost (GeV)

numevents=numevents+(1.*rweight);

hist1nt->Fill(aEvent.ntracks_reco,rweight);
if(aEvent.ccnc_truth==1)
hist2nt->Fill(aEvent.ntracks_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
hist3nt->Fill(aEvent.ntracks_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
{

histeff_qe2->Fill(aEvent.lep_mom_truth,rweight);

hist4nt->Fill(aEvent.ntracks_reco,rweight);
}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
{
hist5nt->Fill(aEvent.ntracks_reco,rweight);
histeff_res2->Fill(aEvent.lep_mom_truth,rweight);
}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
{
hist6nt->Fill(aEvent.ntracks_reco,rweight);
histeff_dis2->Fill(aEvent.lep_mom_truth,rweight);
}

hist1nvt->Fill(aEvent.nvertextracks_reco,rweight);
if(aEvent.ccnc_truth==1)
hist2nvt->Fill(aEvent.nvertextracks_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
hist3nvt->Fill(aEvent.nvertextracks_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
{
hist4nvt->Fill(aEvent.nvertextracks_reco,rweight);
}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
hist5nvt->Fill(aEvent.nvertextracks_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
hist6nvt->Fill(aEvent.nvertextracks_reco,rweight);


hist1ntexit->Fill(aEvent.ntrackendonboundary_reco,rweight);
if(aEvent.ccnc_truth==1)
hist2ntexit->Fill(aEvent.ntrackendonboundary_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
hist3ntexit->Fill(aEvent.ntrackendonboundary_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
{
hist4ntexit->Fill(aEvent.ntrackendonboundary_reco,rweight);
}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
hist5ntexit->Fill(aEvent.ntrackendonboundary_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
hist6ntexit->Fill(aEvent.ntrackendonboundary_reco,rweight);



//not obvious cut
// if(aEvent.nvertextracks_reco<4)
// {

hist1nvtxclusu->Fill(aEvent.nvertexclustersu_reco,rweight);
if(aEvent.ccnc_truth==1)
hist2nvtxclusu->Fill(aEvent.nvertexclustersu_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
hist3nvtxclusu->Fill(aEvent.nvertexclustersu_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
hist4nvtxclusu->Fill(aEvent.nvertexclustersu_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
hist5nvtxclusu->Fill(aEvent.nvertexclustersu_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
hist6nvtxclusu->Fill(aEvent.nvertexclustersu_reco,rweight);

hist1nvtxclusv->Fill(aEvent.nvertexclustersv_reco,rweight);
if(aEvent.ccnc_truth==1)
hist2nvtxclusv->Fill(aEvent.nvertexclustersv_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
hist3nvtxclusv->Fill(aEvent.nvertexclustersv_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
hist4nvtxclusv->Fill(aEvent.nvertexclustersv_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
hist5nvtxclusv->Fill(aEvent.nvertexclustersv_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
hist6nvtxclusv->Fill(aEvent.nvertexclustersv_reco,rweight);



hist1nclusu->Fill(aEvent.nclusu_reco,rweight);
if(aEvent.ccnc_truth==1)
hist2nclusu->Fill(aEvent.nclusu_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
hist3nclusu->Fill(aEvent.nclusu_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
hist4nclusu->Fill(aEvent.nclusu_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
hist5nclusu->Fill(aEvent.nclusu_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
hist6nclusu->Fill(aEvent.nclusu_reco,rweight);

hist1nclusv->Fill(aEvent.nclusv_reco,rweight);
if(aEvent.ccnc_truth==1)
hist2nclusv->Fill(aEvent.nclusv_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
hist3nclusv->Fill(aEvent.nclusv_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
hist4nclusv->Fill(aEvent.nclusv_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
hist5nclusv->Fill(aEvent.nclusv_reco,rweight);
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
hist6nclusv->Fill(aEvent.nclusv_reco,rweight);

//not obvious cut
// if(aEvent.nvertextracks_reco<4) if(aEvent.nvtxclusu_reco<4&&aEvent.nvtxclusv_reco<4) //58, 46
// if(aEvent.nvertexclustersu_reco<3&&aEvent.nvertexclustersv_reco<3
// )
// {


hist11->Fill(aEvent.enu_truth,rweight);
histpur1->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);

if(aEvent.ccnc_truth==1)
hist22->Fill(aEvent.enu_truth,rweight);
if(aEvent.ccnc_truth==0 && aEvent.nuPDG_truth==14)
{
hist33->Fill(aEvent.enu_truth,rweight);
histeff2->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
histeff22->Fill(aEvent.lep_mom_truth,rweight);
histeff222->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
histeff2222->Fill(aEvent.lep_mom_truth,rweight);




histunfold->Fill((180./3.14159)*TMath::ACos(dcosz), (180./3.14159)*TMath::ACos(dcosztruth),rweight);
histunfold2->Fill(((180./3.14159)*TMath::ACos(dcosz)-(180./3.14159)*TMath::ACos(dcosztruth)),rweight);



histunfoldmom->Fill(aEvent.trk_mom_minos+MOM_OFFSET,aEvent.lep_mom_truth,rweight);
histunfoldmom2->Fill(aEvent.trk_mom_minos+MOM_OFFSET-aEvent.lep_mom_truth,rweight);
histunfoldmom22->Fill(aEvent.trk_mom_minos-aEvent.lep_mom_truth,rweight);




int bin=((180./3.14159)*TMath::ACos(dcosz)-0.)*(36./72.)-.0000001;
if(bin<30&&bin>=0)
{


if(((180./3.14159)*TMath::ACos(dcosz)+((180./3.14159)*TMath::ACos(dcosz)*thetasyst[bin]))>0.)
hist3unfoldup->Fill((180./3.14159)*TMath::ACos(dcosz)+((180./3.14159)*TMath::ACos(dcosz)*thetasyst[bin]),(180./3.14159)*TMath::ACos(dcosztruth),rweight);
else
{
hist3unfoldup->Fill(.000001,(180./3.14159)*TMath::ACos(dcosztruth),rweight);
}

hist3unfolddown->Fill((180./3.14159)*TMath::ACos(dcosz)-((180./3.14159)*TMath::ACos(dcosz)*thetasyst[bin]),(180./3.14159)*TMath::ACos(dcosztruth),rweight);



}

int bin2=((aEvent.trk_mom_minos+MOM_OFFSET)-0.)*(40./50.);
if(bin2<30&&bin2>=0)
{
hist3unfoldmueup->Fill((aEvent.trk_mom_minos+MOM_OFFSET)+((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]),aEvent.lep_mom_truth,rweight);

if(((aEvent.trk_mom_minos+MOM_OFFSET)-((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]))>0.)
hist3unfoldmuedown->Fill((aEvent.trk_mom_minos+MOM_OFFSET)-((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]),aEvent.lep_mom_truth,rweight);
else
hist3unfoldmuedown->Fill(.000001,aEvent.lep_mom_truth,rweight);
}











histmomsyst->Fill(aEvent.lep_mom_truth, fabs(aEvent.trk_mom_minos+MOM_OFFSET-aEvent.lep_mom_truth)/aEvent.lep_mom_truth,rweight);

histmomsyst2->Fill(aEvent.lep_mom_truth, (aEvent.trk_mom_minos+MOM_OFFSET-aEvent.lep_mom_truth)/aEvent.lep_mom_truth,rweight);

histmomsyst3->Fill(aEvent.lep_mom_truth, (aEvent.trk_mom_minos+MOM_OFFSET-aEvent.lep_mom_truth)/aEvent.lep_mom_truth,rweight);

histmomsyst2ext->Fill(aEvent.lep_mom_truth, (aEvent.trk_mom_minos+MOM_OFFSET-aEvent.lep_mom_truth)/aEvent.lep_mom_truth,rweight);

histthetasyst->Fill((180./3.14159)*TMath::ACos(dcosztruth),((180./3.14159)*TMath::ACos(dcosz))-((180./3.14159)*TMath::ACos(dcosztruth)),rweight);

histthetasyst2->Fill((180./3.14159)*TMath::ACos(dcosztruth),(((180./3.14159)*TMath::ACos(dcosz))-((180./3.14159)*TMath::ACos(dcosztruth)))/((180./3.14159)*TMath::ACos(dcosztruth)),rweight);

histthetasyst3->Fill((180./3.14159)*TMath::ACos(dcosztruth),(((180./3.14159)*TMath::ACos(dcosz))-((180./3.14159)*TMath::ACos(dcosztruth)))/((180./3.14159)*TMath::ACos(dcosztruth)),rweight);

histthetasyst2ext->Fill((180./3.14159)*TMath::ACos(dcosztruth),(((180./3.14159)*TMath::ACos(dcosz))-((180./3.14159)*TMath::ACos(dcosztruth)))/((180./3.14159)*TMath::ACos(dcosztruth)),rweight);

histpur2->Fill(dcosz,rweight);

histunfold3->Fill((180./3.14159)*TMath::ACos(aEvent.trackstart_dcosx_reco),(180./3.14159)*TMath::ACos(aEvent.lep_dcosx_truth),rweight);
histunfold4->Fill((180./3.14159)*TMath::ACos(aEvent.trackstart_dcosx_reco)-(180./3.14159)*TMath::ACos(aEvent.lep_dcosx_truth),rweight);
histunfold5->Fill((180./3.14159)*TMath::ACos(aEvent.trackstart_dcosy_reco),(180./3.14159)*TMath::ACos(aEvent.lep_dcosy_truth),rweight);
histunfold6->Fill((180./3.14159)*TMath::ACos(aEvent.trackstart_dcosy_reco)-(180./3.14159)*TMath::ACos(aEvent.lep_dcosy_truth),rweight);
histunfold33->Fill((180./3.14159)*TMath::ACos(aEvent.trk_dcosx_minos),(180./3.14159)*TMath::ACos(aEvent.lep_dcosx_truth),rweight);

// std::cout<<TMath::ACos(aEvent.trk_dcosx_minos)<<std::endl;

histunfold44->Fill((180./3.14159)*TMath::ACos(aEvent.trk_dcosx_minos)-(180./3.14159)*TMath::ACos(aEvent.lep_dcosx_truth),rweight);
histunfold55->Fill((180./3.14159)*TMath::ACos(aEvent.trk_dcosy_minos),(180./3.14159)*TMath::ACos(aEvent.lep_dcosy_truth),rweight);
histunfold66->Fill((180./3.14159)*TMath::ACos(aEvent.trk_dcosy_minos)-(180./3.14159)*TMath::ACos(aEvent.lep_dcosy_truth),rweight);





// histtest->Fill(aEvent.trackstart_dcosx_reco,rweight);
// histtest2->Fill(aEvent.trackstart_dcosy_reco,rweight);
// 
// histtest3->Fill(TMath::ATan2(aEvent.trackstart_dcosy_reco,aEvent.trackstart_dcosx_reco),rweight);
// histtest4->Fill(aEvent.trackstart_dcosx_reco, aEvent.trackstart_dcosy_reco,rweight);


vtxhist->Fill(aEvent.vtxx_reco-aEvent.nuvtxx_truth,rweight);
vtxhist2->Fill(aEvent.vtxy_reco-aEvent.nuvtxy_truth,rweight);
vtxhist3->Fill(aEvent.vtxz_reco-aEvent.nuvtxz_truth,rweight);
vtxhistt->Fill(aEvent.vtxx_reco,aEvent.nuvtxx_truth,rweight);
vtxhistt2->Fill(aEvent.vtxy_reco,aEvent.nuvtxy_truth,rweight);
vtxhistt3->Fill(aEvent.vtxz_reco,aEvent.nuvtxz_truth,rweight);

vtxhist4->Fill(sqrt(pow(aEvent.vtxx_reco-aEvent.nuvtxx_truth,2)+pow(aEvent.vtxy_reco-aEvent.nuvtxy_truth,2)+pow(aEvent.vtxz_reco-aEvent.nuvtxz_truth,2)),rweight);

numccevents_aftercuts=numccevents_aftercuts+(1.*rweight);




if(aEvent.lep_mom_truth<25.)
numccevents_aftercuts_mom=numccevents_aftercuts_mom+(1.*rweight);

if((180./3.14159)*TMath::ACos(dcosztruth)<36.)
numccevents_aftercuts_theta=numccevents_aftercuts_theta+(1.*rweight);

}
else
{
histbkgd->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
histbkgd2->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);

}

if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
{
hist44->Fill(aEvent.enu_truth,rweight);
numccqeevents_aftercuts=numccqeevents_aftercuts+(1.*rweight);
}

if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
{
hist55->Fill(aEvent.enu_truth,rweight);
numccresevents_aftercuts=numccresevents_aftercuts+(1.*rweight);
}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
{
hist66->Fill(aEvent.enu_truth,rweight);
numccdisevents_aftercuts=numccdisevents_aftercuts+(1.*rweight);
}


// }



 // }


}

// hist1->Fill((180./3.14159)*TMath::ACos(aEvent.trackstart_dcosz_reco), (180./3.14159)*TMath::ACos(aEvent.lep_dcosz_truth));
// 
// hist2->Fill((180./3.14159)*TMath::ACos(aEvent.lep_dcosz_reco));
// hist3->Fill((180./3.14159)*TMath::ACos(aEvent.trackstart_dcosz_reco));

// std::cout<<aEvent.trackstart_dcosz_reco<<" "<<aEvent.lep_dcosz_truth<<std::endl;

}

}
std::cout<<numccevents_beforecuts<<" "<<numccevents_aftercuts<<std::endl;

double ccminosefficiency_beforecuts=(double)numccevents_beforecuts_recminos/(double)numccevents_beforecuts;

double ccqeminosefficiency_beforecuts=(double)numccqeevents_beforecuts_recminos/(double)numccqeevents_beforecuts;

double ccresminosefficiency_beforecuts=(double)numccresevents_beforecuts_recminos/(double)numccresevents_beforecuts;

double ccdisminosefficiency_beforecuts=(double)numccdisevents_beforecuts_recminos/(double)numccdisevents_beforecuts;

double ccminosefficiency_aftercuts=(double)numccevents_aftercuts_recminos/(double)numccevents_beforecuts;



double ccefficiency=(double)numccevents_aftercuts/(double)numccevents_beforecuts;

double ccefficiency_theta=(double)numccevents_aftercuts_theta/(double)numccevents_beforecuts_theta;
double ccefficiency_mom=(double)numccevents_aftercuts_mom/(double)numccevents_beforecuts_mom;

double ccpurity=(double)numccevents_aftercuts/(double)numevents;
std::cout<<"CC Efficiency "<<ccefficiency<<" "<<numccevents_aftercuts<<" "<<numccevents_beforecuts<<std::endl;

std::cout<<"CC Efficiency theta "<<ccefficiency_theta<<" "<<numccevents_aftercuts_theta<<" "<<numccevents_beforecuts_theta<<std::endl;

std::cout<<"CC Efficiency mom "<<ccefficiency_mom<<" "<<numccevents_aftercuts_mom<<" "<<numccevents_beforecuts_mom<<std::endl;

double ccqeefficiency=(double)numccqeevents_aftercuts/(double)numccqeevents_beforecuts;
std::cout<<"CCQE Efficiency "<<ccqeefficiency<<std::endl;

double ccresefficiency=(double)numccresevents_aftercuts/(double)numccresevents_beforecuts;
std::cout<<"CCRES Efficiency "<<ccresefficiency<<std::endl;

double ccdisefficiency=(double)numccdisevents_aftercuts/(double)numccdisevents_beforecuts;
std::cout<<"CCDIS Efficiency "<<ccdisefficiency<<std::endl;


std::cout<<"CC MINOS Efficiency before cuts "<<ccminosefficiency_beforecuts<<std::endl;

std::cout<<"CCQE MINOS Efficiency before cuts "<<ccqeminosefficiency_beforecuts<<std::endl;

std::cout<<"CCRES MINOS Efficiency before cuts "<<ccresminosefficiency_beforecuts<<std::endl;

std::cout<<"CCDIS MINOS Efficiency before cuts "<<ccdisminosefficiency_beforecuts<<std::endl;

std::cout<<"CC MINOS Efficiency after cuts "<<ccminosefficiency_aftercuts<<std::endl;




std::cout<<"ArgoNeuT+Matching CC efficiency "<<ccefficiency/ccminosefficiency_beforecuts<<std::endl;


std::cout<<"ArgoNeuT+Matching CCQE efficiency "<<ccqeefficiency/ccqeminosefficiency_beforecuts<<std::endl;

std::cout<<"ArgoNeuT+Matching CCRES efficiency "<<ccresefficiency/ccresminosefficiency_beforecuts<<std::endl;

std::cout<<"ArgoNeuT+Matching CCDIS efficiency "<<ccdisefficiency/ccdisminosefficiency_beforecuts<<std::endl;

std::cout<<numncevents<<" "<<numtotalevents<<std::endl;

// new TCanvas;
// hist1->SetMarkerStyle(34);
// hist1->SetTitle("Neutrino candidate vertex positions (neutrino-mode)");
// hist1->SetXTitle("Z coordinate (cm)");
// hist1->SetYTitle("X coordinate (cm)");
// hist1->DrawCopy("colz");

char exxpression[64], exxpressionn[64];
sprintf(exxpression,"After correction");
TLatex* laex = new TLatex(.2, .8, exxpression);
laex->SetNDC();
laex->SetTextSize(.038);
laex->SetTextColor(1);

sprintf(exxpressionn,"Before correction");
TLatex* laexx = new TLatex(.2, .75, exxpressionn);
laexx->SetNDC();
laexx->SetTextSize(.038);
laexx->SetTextColor(2);


char expression[64], expressionn[64], expressionnn[64], expressionnnn[64], expressionnnnn[64], expressionnnnnn[64];
sprintf(expression,"All flavors/interactions");
TLatex* latex = new TLatex(.5, .8, expression);
latex->SetNDC();
latex->SetTextSize(.038);
latex->SetTextColor(1);

sprintf(expressionn,"NC/WS");
TLatex* latexx = new TLatex(.5, .75, expressionn);
latexx->SetNDC();
latexx->SetTextSize(.038);
latexx->SetTextColor(2);


sprintf(expressionnn,"#nu_{#mu} CC");
TLatex* latexxx = new TLatex(.5, .7, expressionnn);
latexxx->SetNDC();
latexxx->SetTextSize(.038);
latexxx->SetTextColor(3);

sprintf(expressionnnn,"#nu_{#mu} CCQE");
TLatex* latexxxx = new TLatex(.5, .65, expressionnnn);
latexxxx->SetNDC();
latexxxx->SetTextSize(.038);
latexxxx->SetTextColor(4);

sprintf(expressionnnnn,"#nu_{#mu} CCRES");
TLatex* latexxxxx = new TLatex(.5, .6, expressionnnnn);
latexxxxx->SetNDC();
latexxxxx->SetTextSize(.038);
latexxxxx->SetTextColor(7);

sprintf(expressionnnnnn,"#nu_{#mu} CCDIS");
TLatex* latexxxxxx = new TLatex(.5, .55, expressionnnnnn);
latexxxxxx->SetNDC();
latexxxxxx->SetTextSize(.038);
latexxxxxx->SetTextColor(6);


double potfactor=8.518840e18/totpot;

TLine l;



new TCanvas;


gStyle->SetOptStat(0);

 hist1->SetTitle("Neutrino energy truth (before cuts)");
 hist1->SetXTitle("E_{#nu} (GeV)");
 hist1->SetYTitle("events");
 hist1->DrawNormalized("", hist1->Integral()*potfactor);
 hist2->SetLineColor(2);
 hist2->DrawNormalized("SAME", hist2->Integral()*potfactor);
 hist3->SetLineColor(3);
 hist3->DrawNormalized("SAME", hist3->Integral()*potfactor);
 hist4->SetLineColor(4);
 hist4->DrawNormalized("SAME", hist4->Integral()*potfactor);
 hist5->SetLineColor(7);
 hist5->DrawNormalized("SAME", hist5->Integral()*potfactor);
 hist6->SetLineColor(6);
 hist6->DrawNormalized("SAME", hist6->Integral()*potfactor);
latex->Draw();
latexx->Draw();
latexxx->Draw();
latexxxx->Draw();
latexxxxx->Draw();
latexxxxxx->Draw();
 gPad->GetCanvas()->Print("histnuE_beforecuts.pdf"); 
 


 
 
new TCanvas;
 hist11->SetTitle("Neutrino energy truth (after cuts)");
 hist11->SetXTitle("E_{#nu} (GeV)");
 hist11->SetYTitle("events");
 hist11->DrawNormalized("", hist11->Integral()*potfactor);
 hist22->SetLineColor(2);
 hist22->DrawNormalized("SAME", hist22->Integral()*potfactor);
 hist33->SetLineColor(3);
 hist33->DrawNormalized("SAME", hist33->Integral()*potfactor);
 hist44->SetLineColor(4);
 hist44->DrawNormalized("SAME", hist44->Integral()*potfactor);
 hist55->SetLineColor(7);
 hist55->DrawNormalized("SAME", hist55->Integral()*potfactor);
 hist66->SetLineColor(6);
 hist6->DrawNormalized("SAME", hist66->Integral()*potfactor);
 latex->Draw();
latexx->Draw();
latexxx->Draw();
latexxxx->Draw();
latexxxxx->Draw();
latexxxxxx->Draw();
 gPad->GetCanvas()->Print("histnuE_aftercuts.pdf");
 
 
 
 
new TCanvas;
gStyle->SetOptStat();
 hist1nt->SetTitle("Number of 3D tracks (after cuts)");
 hist1nt->SetXTitle("number of 3D tracks");
 hist1nt->SetYTitle("events");
 hist1nt->DrawNormalized("", hist1nt->Integral()*potfactor);
 hist2nt->SetLineColor(2);
 hist2nt->DrawNormalized("SAME", hist2nt->Integral()*potfactor);
 hist3nt->SetLineColor(3);
 hist3nt->DrawNormalized("SAME", hist3nt->Integral()*potfactor);
 hist4nt->SetLineColor(4);
 hist4nt->DrawNormalized("SAME", hist4nt->Integral()*potfactor);
 hist5nt->SetLineColor(7);
 hist5nt->DrawNormalized("SAME", hist5nt->Integral()*potfactor);
 hist6nt->SetLineColor(6);
 hist6nt->DrawNormalized("SAME", hist6nt->Integral()*potfactor);
 
 latex->Draw();
latexx->Draw();
latexxx->Draw();
latexxxx->Draw();
latexxxxx->Draw();
latexxxxxx->Draw();

// new TCanvas;
//  hist1nvt->SetTitle("Number of 3D tracks associated with vertex (after cuts)");
//  hist1nvt->SetXTitle("number of 3D tracks");
//  hist1nvt->SetYTitle("events");
//  hist1nvt->DrawNormalized("", hist1nvt->Integral()*potfactor);
//  hist2nvt->SetLineColor(2);
//  hist2nvt->DrawNormalized("SAME",hist2nvt->Integral()*potfactor);
//  hist3nvt->SetLineColor(3);
//  hist3nvt->DrawNormalized("SAME",hist3nvt->Integral()*potfactor);
//  hist4nvt->SetLineColor(4);
//  hist4nvt->DrawNormalized("SAME",hist4nvt->Integral()*potfactor);
//  hist5nvt->SetLineColor(7);
//  hist5nvt->DrawNormalized("SAME",hist5nvt->Integral()*potfactor);
//  hist6nvt->SetLineColor(6);
//  hist6nvt->DrawNormalized("SAME",hist6nvt->Integral()*potfactor);
//  
//  latex->Draw();
// latexx->Draw();
// latexxx->Draw();
// latexxxx->Draw();
// latexxxxx->Draw();
// latexxxxxx->Draw();

new TCanvas;
 hist1ntexit->SetTitle("Number of exiting 3D tracks (after cuts)");
 hist1ntexit->SetXTitle("number of exiting 3D tracks");
 hist1ntexit->SetYTitle("events");
 hist1ntexit->DrawNormalized("", hist1ntexit->Integral()*potfactor );
 hist2ntexit->SetLineColor(2);
 hist2ntexit->DrawNormalized("SAME", hist2ntexit->Integral()*potfactor);
 hist3ntexit->SetLineColor(3);
 hist3ntexit->DrawNormalized("SAME", hist3ntexit->Integral()*potfactor);
 hist4ntexit->SetLineColor(4);
 hist4ntexit->DrawNormalized("SAME", hist4ntexit->Integral()*potfactor);
 hist5ntexit->SetLineColor(7);
 hist5ntexit->DrawNormalized("SAME", hist5ntexit->Integral()*potfactor);
 hist6ntexit->SetLineColor(6);
 hist6ntexit->DrawNormalized("SAME", hist6ntexit->Integral()*potfactor);

 latex->Draw();
latexx->Draw();
latexxx->Draw();
latexxxx->Draw();
latexxxxx->Draw();
latexxxxxx->Draw();

 new TCanvas;
 hist1nvtxclusu->SetTitle("Number of 2D merged-lines associated with vertex in collection plane (after cuts)");
 hist1nvtxclusu->SetXTitle("number of merged-lines");
 hist1nvtxclusu->SetYTitle("events");
 hist1nvtxclusu->DrawNormalized("",hist1nvtxclusu->Integral()*potfactor);
 hist2nvtxclusu->SetLineColor(2);
 hist2nvtxclusu->DrawNormalized("SAME",hist2nvtxclusu->Integral()*potfactor);
 hist3nvtxclusu->SetLineColor(3);
 hist3nvtxclusu->DrawNormalized("SAME",hist3nvtxclusu->Integral()*potfactor);
 hist4nvtxclusu->SetLineColor(4);
 hist4nvtxclusu->DrawNormalized("SAME",hist4nvtxclusu->Integral()*potfactor);
 hist5nvtxclusu->SetLineColor(7);
 hist5nvtxclusu->DrawNormalized("SAME",hist5nvtxclusu->Integral()*potfactor);
 hist6nvtxclusu->SetLineColor(6);
 hist6nvtxclusu->DrawNormalized("SAME",hist6nvtxclusu->Integral()*potfactor);

 latex->Draw();
latexx->Draw();
latexxx->Draw();
latexxxx->Draw();
latexxxxx->Draw();
latexxxxxx->Draw();

 new TCanvas;
 hist1nvtxclusv->SetTitle("Number of 2D merged-lines associated with vertex in induction plane (after cuts)");
 hist1nvtxclusv->SetXTitle("number of merged-lines");
 hist1nvtxclusv->SetYTitle("events");
 hist1nvtxclusv->DrawNormalized("",hist1nvtxclusv->Integral()*potfactor);
 hist2nvtxclusv->SetLineColor(2);
 hist2nvtxclusv->DrawNormalized("SAME",hist2nvtxclusv->Integral()*potfactor);
 hist3nvtxclusv->SetLineColor(3);
 hist3nvtxclusv->DrawNormalized("SAME",hist3nvtxclusv->Integral()*potfactor);
 hist4nvtxclusv->SetLineColor(4);
 hist4nvtxclusv->DrawNormalized("SAME",hist4nvtxclusv->Integral()*potfactor);
 hist5nvtxclusv->SetLineColor(7);
 hist5nvtxclusv->DrawNormalized("SAME",hist5nvtxclusv->Integral()*potfactor);
 hist6nvtxclusv->SetLineColor(6);
 hist6nvtxclusv->DrawNormalized("SAME",hist6nvtxclusv->Integral()*potfactor);
  latex->Draw();
latexx->Draw();
latexxx->Draw();
latexxxx->Draw();
latexxxxx->Draw();
latexxxxxx->Draw();
 
 
  new TCanvas;
 hist1nclusu->SetTitle("Number of 2D merged-lines in collection plane (after cuts)");
 hist1nclusu->SetXTitle("number of merged-lines");
 hist1nclusu->SetYTitle("events");
 hist1nclusu->DrawNormalized("",hist1nclusu->Integral()*potfactor);
 hist2nclusu->SetLineColor(2);
 hist2nclusu->DrawNormalized("SAME",hist2nclusu->Integral()*potfactor);
 hist3nclusu->SetLineColor(3);
 hist3nclusu->DrawNormalized("SAME",hist3nclusu->Integral()*potfactor);
 hist4nclusu->SetLineColor(4);
 hist4nclusu->DrawNormalized("SAME",hist4nclusu->Integral()*potfactor);
 hist5nclusu->SetLineColor(7);
 hist5nclusu->DrawNormalized("SAME",hist5nclusu->Integral()*potfactor);
 hist6nclusu->SetLineColor(6);
 hist6nclusu->DrawNormalized("SAME",hist6nclusu->Integral()*potfactor);
 
 latex->Draw();
latexx->Draw();
latexxx->Draw();
latexxxx->Draw();
latexxxxx->Draw();
latexxxxxx->Draw();

 new TCanvas;
 hist1nclusv->SetTitle("Number of 2D merged-lines in induction plane (after cuts)");
 hist1nclusv->SetXTitle("number of merged-lines");
 hist1nclusv->SetYTitle("events");
 hist1nclusv->DrawNormalized("",hist1nclusv->Integral()*potfactor);
 hist2nclusv->SetLineColor(2);
 hist2nclusv->DrawNormalized("SAME",hist2nclusv->Integral()*potfactor);
 hist3nclusv->SetLineColor(3);
 hist3nclusv->DrawNormalized("SAME",hist3nclusv->Integral()*potfactor);
 hist4nclusv->SetLineColor(4);
 hist4nclusv->DrawNormalized("SAME",hist4nclusv->Integral()*potfactor);
 hist5nclusv->SetLineColor(7);
 hist5nclusv->DrawNormalized("SAME",hist5nclusv->Integral()*potfactor);
 hist6nclusv->SetLineColor(6);
 hist6nclusv->DrawNormalized("SAME",hist6nclusv->Integral()*potfactor);
 
 latex->Draw();
latexx->Draw();
latexxx->Draw();
latexxxx->Draw();
latexxxxx->Draw();
latexxxxxx->Draw();

// 
// l.DrawLine(0,0,50,50);
// 



 new TCanvas;
 
 gStyle->SetOptStat(0);
 
 histeff2->Divide(histeff1);
 histeff2->SetTitle("CC #nu_{#mu} muon reconstruction probability (after cuts)");
 histeff2->SetXTitle("#theta_{#mu, true} (degrees) ");
 histeff2->SetYTitle("probability");
 histeff2->SetMaximum(1);
 histeff2->SetMinimum(0);
 histeff2->DrawCopy("E1");
gPad->GetCanvas()->Print("histeff2.pdf");

 for(int i=1;i<= histeff2->GetNbinsX();i++)
{
std::cout<<histeff2->GetBinCenter(i)<<" "<<histeff2->GetBinContent(i)<<std::endl;

}

 new TCanvas;
 histeff22->Divide(histeff11);
 histeff22->SetTitle("CC #nu_{#mu} muon reconstruction probability (after cuts)");
 histeff22->SetXTitle("P_{#mu, true} (GeV/c) ");
 histeff22->SetYTitle("probability");
 histeff22->SetMaximum(1);
 histeff22->SetMinimum(0);
 histeff22->DrawCopy("E1");
gPad->GetCanvas()->Print("histeff22.pdf");

 new TCanvas;
//  histeff2->Multiply(histeff1);
  histeff222->Divide(histeff111);
 histeff222->SetTitle("CC #nu_{#mu} muon ArgoNeuT+matching reconstruction probability (after cuts)");
 histeff222->SetXTitle("#theta_{#mu, true} (degrees)");
 histeff222->SetYTitle("probability");
 histeff222->SetMaximum(1);
 histeff222->SetMinimum(0);
 histeff222->DrawCopy("E1");
gPad->GetCanvas()->Print("histeff222.pdf");


 new TCanvas;
 histeff2222->Divide(histeff1111);
 histeff2222->SetTitle("CC #nu_{#mu} muon ArgoNeuT+matching reconstruction probability (after cuts)");
 histeff2222->SetXTitle("P_{#mu, true} (GeV/c) ");
 histeff2222->SetYTitle("probability");
 histeff2222->SetMaximum(1);
 histeff2222->SetMinimum(0);
 histeff2222->DrawCopy("E1");
 gPad->GetCanvas()->Print("histeff2222.pdf");
 
 
 new TCanvas;
 histeff111->Divide(histeff1);
 histeff111->SetTitle("CC #nu_{#mu} muon MINOS reconstruction probability (after cuts)");
 histeff111->SetXTitle("#theta_{#mu, true}");
 histeff111->SetYTitle("probability");
 histeff111->SetMaximum(1);
 histeff111->SetMinimum(0);
 histeff111->Draw();
   gPad->GetCanvas()->Print("histminoseff.pdf");
  new TCanvas;
 histeff1111->Divide(histeff11);
 histeff1111->SetTitle("CC #nu_{#mu} muon MINOS reconstruction probability (after cuts)");
 histeff1111->SetXTitle("P_{#mu} (GeV/c)");
 histeff1111->SetYTitle("probability");
 histeff1111->SetMaximum(1);
 histeff1111->SetMinimum(0);
 histeff1111->Draw();
   gPad->GetCanvas()->Print("histminoseff2.pdf");
 new TCanvas;
 histpur2->Divide(histpur1);
 histpur2->SetTitle("CC #nu_{#mu} muon reconstruction purity (after cuts)");
 histpur2->SetXTitle("#theta_{#mu, recosim} ");
 histpur2->SetYTitle("purity after cuts & matching");
 histpur2->DrawCopy("E1");
 

 
 new TCanvas;
 histunfold->SetTitle("CC #nu_{#mu} muon #theta_{#mu} recosim and truth (after cuts)");
 histunfold->SetXTitle("#theta_{#mu, recosim} (degrees)");
 histunfold->SetYTitle("#theta_{#mu, true} (degrees)");
 histunfold->SetMinimum(histunfold->GetMaximum()*.0);
 histunfold->GetYaxis()->SetTitleOffset(1.3);
 gPad->GetCanvas()->SetLogz();
 histunfold->GetXaxis()->SetRangeUser(0.,35.999);
 histunfold->GetYaxis()->SetRangeUser(0.,35.999);
 histunfold->DrawCopy("colz");
 histunfold->GetXaxis()->SetRangeUser(0.,35.999);
 histunfold->GetYaxis()->SetRangeUser(0.,35.999);
 
 
 
 TFile *unf=new TFile("unfold.root", "RECREATE"); 
 histunfold->Write();
 hist3unfoldup->Write();
 hist3unfolddown->Write();
 unf->Close();
 
 TFile *en=new TFile("en.root", "RECREATE"); 
 histunweighted->Scale(potfactor);
 histunweighted->Write();
 en->Close();
 
 l.DrawLine(0,0,36,36);
  gPad->GetCanvas()->Print("histunfold.pdf");
 new TCanvas;
 histunfold->SetTitle("CC #nu_{#mu} muon #theta_{#mu} recosim and truth (after cuts)");
 histunfold->SetXTitle("#theta_{#mu,recosim} (degrees)");
 histunfold->SetYTitle("#theta_{#mu, true} (degrees)");
//  histunfold->GetXaxis()->SetRangeUser(.8,1.);
//  histunfold->GetYaxis()->SetRangeUser(.8,1.);
 histunfold->GetXaxis()->SetTitleOffset(1.9);
 histunfold->GetYaxis()->SetTitleOffset(1.9);
 histunfold->GetXaxis()->SetRangeUser(0.,36.);
 histunfold->GetYaxis()->SetRangeUser(0.,36.);
 histunfold->DrawCopy("lego2");
  histunfold->GetXaxis()->SetRangeUser(0.,36.);
 histunfold->GetYaxis()->SetRangeUser(0.,36.);
//  l.DrawLine(.8,.8,1,1);
 gPad->GetCanvas()->Print("histunfoldzoomed.pdf");
 
 
  new TCanvas;

 

 histunfold2->SetTitle("CC #nu_{#mu} muon #theta_{#mu} recosim - truth (after cuts)");
 histunfold2->SetXTitle("#theta_{#mu,recosim} - #theta_{#mu,true} (degrees)");
 histunfold2->SetYTitle("fraction of events");
 histunfold2->GetYaxis()->SetTitleOffset(1.27);
 histunfold2->Scale(1./histunfold2->Integral());
 histunfold2->Fit("gaus","W","",-2.5,2.5);
  gStyle->SetOptStat(1000001110);
gPad->GetCanvas()->Print("histunfold2.pdf");
 
  new TCanvas;
   gStyle->SetOptStat(0);
 histunfoldmom->SetTitle("CC #nu_{#mu} muon momentum recosim and truth (after cuts)");
 histunfoldmom->SetXTitle("P_{#mu, recosim} (GeV/c)");
 histunfoldmom->SetYTitle("P_{#mu, true} (GeV/c)");
 histunfoldmom->SetMinimum(histunfoldmom->GetMaximum()*.0);
 gPad->GetCanvas()->SetLogz();
  histunfoldmom->GetXaxis()->SetRangeUser(0.,24.999);
 histunfoldmom->GetYaxis()->SetRangeUser(0.,24.999);
 histunfoldmom->DrawCopy("colz");
  histunfoldmom->GetXaxis()->SetRangeUser(0.,24.999);
 histunfoldmom->GetYaxis()->SetRangeUser(0.,24.999);
 
 
 
  TFile *unfmom=new TFile("unfoldmom.root", "RECREATE"); 
 histunfoldmom->Write();
 hist3unfoldmueup->Write();
 hist3unfoldmuedown->Write(); 
 unfmom->Close();
 
 
 l.DrawLine(0,0,25,25);
 gPad->GetCanvas()->Print("histunfoldmom.pdf");

  new TCanvas;
   gStyle->SetOptStat(0);
 histunfoldmom->SetTitle("CC #nu_{#mu} muon momentum recosim and truth (after cuts)");
 histunfoldmom->SetXTitle("P_{#mu, recosim} (GeV/c)");
 histunfoldmom->SetYTitle("P_{#mu, true} (GeV/c)");
 histunfoldmom->GetXaxis()->SetRangeUser(0.,24.999);
 histunfoldmom->GetYaxis()->SetRangeUser(0.,24.999);
//  histunfoldmom->SetMinimum(histunfoldmom->GetMaximum()*.000833);
 histunfoldmom->GetXaxis()->SetTitleOffset(1.6);
 histunfoldmom->GetYaxis()->SetTitleOffset(1.6);
 histunfoldmom->DrawCopy("lego2");
   histunfoldmom->GetXaxis()->SetRangeUser(0.,24.999);
 histunfoldmom->GetYaxis()->SetRangeUser(0.,24.999);
//  l.DrawLine(0,0,25,25);
 gPad->GetCanvas()->Print("histunfoldmomzoomed.pdf");


 

 new TCanvas;
  gStyle->SetOptStat(0);
 histunfold3->SetTitle("CC #nu_{#mu} muon #theta_{x} recosim and truth (after cuts)");
 histunfold3->SetXTitle("#theta_{x, recosim} (degrees)");
 histunfold3->SetYTitle("#theta_{x, true} (degrees)");
 histunfold3->SetMinimum(histunfold3->GetMaximum()*.01);
 gPad->GetCanvas()->SetLogz(); 
 histunfold3->DrawCopy("colz");
 histunfold3->GetXaxis()->SetRangeUser(0.,36.);
 histunfold3->GetYaxis()->SetRangeUser(0.,36.);
  l.DrawLine(0,0,1.,1.);
 gPad->GetCanvas()->Print("histunfold3.pdf");
 


 
new TCanvas;

 histunfold4->SetTitle("CC #nu_{#mu} muon #theta_{x} recosim - truth (after cuts)");
 histunfold4->SetXTitle("#theta_{x, recosim} - #theta_{x, true} (degrees)");
 histunfold4->GetYaxis()->SetTitleOffset(1.27);
 histunfold4->SetYTitle("fraction of events");
 histunfold4->DrawNormalized();
    gStyle->SetOptStat(1000001110);
 gPad->GetCanvas()->Print("histunfold4.pdf");


 
 new TCanvas;
  gStyle->SetOptStat(0);
 histunfold5->SetTitle("CC #nu_{#mu} muon #theta_{y} recosim and truth (after cuts)");
 histunfold5->SetXTitle("#theta_{y, recosim} (degrees)");
 histunfold5->SetYTitle("#theta_{y, true} (degrees)");
 histunfold5->SetMinimum(20);
 gPad->GetCanvas()->SetLogz();
 histunfold5->DrawCopy("colz");
  l.DrawLine(0,0,1.,1.);
  
 
  
 gPad->GetCanvas()->Print("histunfold5.pdf");
 new TCanvas;

 histunfold6->SetTitle("CC #nu_{#mu} muon #theta_{y} recosim - truth (after cuts)");
 histunfold6->SetXTitle("#theta_{y, recosim} - #theta_{y, true} (degrees)");
 histunfold6->GetYaxis()->SetTitleOffset(1.27);
 histunfold6->SetYTitle("fraction of events");
 histunfold6->DrawNormalized(); 
    gStyle->SetOptStat(1000001110);
 gPad->GetCanvas()->Print("histunfold6.pdf");


 
 new TCanvas;
  gStyle->SetOptStat(0);
 histunfold33->SetTitle("CC #nu_{#mu} muon #theta_{x} recominos and truth (after cuts)");
 histunfold33->SetXTitle("#theta_{x, recominos} (degrees)");
 histunfold33->SetYTitle("#theta_{x, true} (degrees)");
 histunfold33->SetMinimum(histunfold33->GetMaximum()*.01);
 histunfold33->DrawCopy("colz");
 l.DrawLine(0,0,1.,1.);
 gPad->GetCanvas()->Print("histunfold33.pdf");
 

 
new TCanvas;

 histunfold44->SetTitle("CC #nu_{#mu} muon #theta_{x} recominos - truth (after cuts)");
 histunfold44->SetXTitle("#theta_{x, recominos} - #theta_{x, true} (degrees)");
 histunfold44->GetYaxis()->SetTitleOffset(1.27);
 histunfold44->SetYTitle("fraction of events");
 histunfold44->DrawNormalized();
    gStyle->SetOptStat(1000001110);
 gPad->GetCanvas()->Print("histunfold44.pdf");


 
 new TCanvas;
  gStyle->SetOptStat(0);
 histunfold55->SetTitle("CC #nu_{#mu} muon #theta_{y} recominos and truth (after cuts)");
 histunfold55->SetXTitle("#theta_{y, recominos} (degrees)");
 histunfold55->SetYTitle("#theta_{y, true} (degrees)");
 histunfold55->SetMinimum(histunfold55->GetMaximum()*.01);
 histunfold55->DrawCopy("colz");
 l.DrawLine(0,0,1.,1.);
 gPad->GetCanvas()->Print("histunfold55.pdf");
 

 
new TCanvas;

 histunfold66->SetTitle("CC #nu_{#mu} muon #theta_{y} recominos - truth (after cuts)");
 histunfold66->SetXTitle("#theta_{y, recominos} - #theta_{y, true} (degrees)");
 histunfold66->GetYaxis()->SetTitleOffset(1.27);
 histunfold66->SetYTitle("fraction of events");
 histunfold66->DrawNormalized();
   gStyle->SetOptStat(1000001110);
 gPad->GetCanvas()->Print("histunfold66.pdf");


 new TCanvas;
 vtxhist->SetTitle("CC #nu_{#mu} X vertex recosim - truth (after cuts)");
 vtxhist->SetXTitle("X recosim - X truth (cm)");
 vtxhist->GetYaxis()->SetTitleOffset(1.27);
 vtxhist->SetYTitle("fraction of events"); 
 vtxhist->Scale(1./vtxhist->Integral());
 vtxhist->Fit("gaus","W","",-10,10);
 gPad->GetCanvas()->Print("vtxhist.pdf");
  new TCanvas;
 vtxhist2->SetTitle("CC #nu_{#mu} Y vertex recosim - truth (after cuts)");
 vtxhist2->SetXTitle("Y recosim - Y truth (cm)");
 vtxhist2->SetYTitle("fraction of events");
 vtxhist2->GetYaxis()->SetTitleOffset(1.27);
 vtxhist2->Scale(1./vtxhist2->Integral());
 vtxhist2->Fit("gaus","W","",-10,10);
 gPad->GetCanvas()->Print("vtxhist2.pdf");
  new TCanvas;
 vtxhist3->SetTitle("CC #nu_{#mu} Z vertex recosim - truth (after cuts)");
 vtxhist3->SetXTitle("Z recosim - Z truth (cm)");
 vtxhist3->GetYaxis()->SetTitleOffset(1.27);
 vtxhist3->SetYTitle("fraction of events"); 
 vtxhist3->Scale(1./vtxhist3->Integral());
  vtxhist3->Fit("gaus","W","",-10,10);
 gPad->GetCanvas()->Print("vtxhist3.pdf");
 

 
 new TCanvas;
  gStyle->SetOptStat(0);
 vtxhistt->SetTitle("CC #nu_{#mu} X vertex recosim and truth (after cuts)");
 vtxhistt->SetXTitle("X recosim (cm)");
 vtxhistt->SetYTitle("X truth (cm)");
 vtxhistt->SetMinimum(vtxhistt->GetMaximum()*.01);
 vtxhistt->DrawCopy("colz");
  l.DrawLine(0,0,47,47);
 gPad->GetCanvas()->Print("vtxhistt.pdf"); 
  new TCanvas;
 vtxhistt2->SetTitle("CC #nu_{#mu} Y vertex recosim and truth (after cuts)");
 vtxhistt2->SetXTitle("Y recosim (cm)");
 vtxhistt2->SetYTitle("Y truth (cm)");
 vtxhistt2->SetMinimum(vtxhistt2->GetMaximum()*.01);
 vtxhistt2->DrawCopy("colz");
  l.DrawLine(-20,-20,20,20);
  gPad->GetCanvas()->Print("vtxhistt2.pdf");
  new TCanvas;
 vtxhistt3->SetTitle("CC #nu_{#mu} Z vertex recosim and truth (after cuts)");
 vtxhistt3->SetXTitle("Z recosim (cm)");
 vtxhistt3->SetYTitle("Z truth (cm)"); 
 vtxhistt3->SetMinimum(vtxhistt3->GetMaximum()*.01);
 vtxhistt3->DrawCopy("colz");
  l.DrawLine(0,0,90,90);
gPad->GetCanvas()->Print("vtxhistt3.pdf");




TProfile::Approximate(kTRUE);  
   new TCanvas;
   TProfile *profmomsyst = histmomsyst->ProfileX("profmomsyst", 0, 50, ""); 
   profmomsyst->SetTitle("CC #nu_{#mu} muon momentum fractional error (comparing recosim & truth)");
   profmomsyst->SetXTitle("P_{true}");
   profmomsyst->SetYTitle("|#DeltaP|/P_{true}");
   profmomsyst->GetYaxis()->SetTitleOffset(1.27);
   // profmomsyst->DrawCopy("E1");
//  gPad->GetCanvas()->Print("profmomsyst.pdf");

     new TCanvas;
     histthetasyst->DrawCopy("colz");
    
   new TCanvas;
   TProfile *profthetasyst = histthetasyst->ProfileX("profthetasyst", 0, 50, ""); 
   profthetasyst->SetTitle("CC #nu_{#mu} muon #theta_{#mu} fractional error (comparing recosim & truth)");
   profthetasyst->SetXTitle("#theta_{#mu, true} (degrees)");
   profthetasyst->SetYTitle("#Delta#theta_{#mu}");
   profthetasyst->GetYaxis()->SetTitleOffset(1.27);
    profthetasyst->DrawCopy("E1");
//gPad->GetCanvas()->Print("profthetasyst.pdf");



TH1D *h2 = histthetasyst->ProjectionY("h2", 1, 1);

TH1D* h2final = new TH1D("h2final", "",18,0.,36);

for(int i=1;i<=histthetasyst->GetNbinsX();i++)
{
 //new TCanvas;
h2 = histthetasyst->ProjectionY("h2", i, i); // where firstYbin = 0 and lastYbin = 9

double minf=-5., maxf=5.;
TF1 *gfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
h2->Fit("Gaussian","RQN"); 
Double_t sigma = gfit->GetParameter(2); //value of 1st parameter
Double_t esigma = gfit->GetParError(2); //error on 1st parameter
h2final->Fill(h2final->GetBinCenter(i),sigma);
h2final->SetBinError(i,esigma);
// py->Draw();
}

new TCanvas;
 gStyle->SetOptStat(0);
h2final->SetTitle("CC #nu_{#mu} muon #theta_{#mu} error (comparing recosim & truth)");
h2final->SetYTitle("#delta[#Delta#theta_{#mu}]");
h2final->SetXTitle("#theta_{#mu, true}");
h2final->GetYaxis()->SetTitleOffset(1.27);
h2final->SetMinimum(0);
h2final->Draw();
gStyle->SetOptStat(0);
gPad->GetCanvas()->Print("thetaerror.pdf");


//    new TCanvas; 
//    profthetasyst->GetXaxis()->SetRangeUser(0.,1.);
//    profthetasyst->GetYaxis()->SetLabelSize(.035);
//    profthetasyst->DrawCopy("E1");
//gPad->GetCanvas()->Print("profthetasyst2.pdf");


 gStyle->SetOptStat(0);

 new TCanvas;
 
 TProfile *profthetasystt = histthetasyst2ext->ProfileX("profthetasystt", 0, 100, "");
 

 histbkgd->SetTitle("CC #nu_{#mu} muon background (after cuts)");
 histbkgd->SetXTitle("#theta_{#mu, recosim} (degrees)");
 histbkgd->SetYTitle("background events");
 histbkgd->Scale(potfactor);
 histbkgd->Draw();
 
 new TCanvas;
 histbkgd2->SetTitle("CC #nu_{#mu} muon background (after cuts)");
 histbkgd2->SetXTitle("P_{#mu, recosim} (GeV/c)");
 histbkgd2->SetYTitle("background events");
 histbkgd2->Scale(potfactor);
 histbkgd2->Draw();

new TCanvas;
histdistance->Draw();// 
  gStyle->SetOptStat(0);
new TCanvas;
  gStyle->SetOptStat(0);
histmomsyst3->SetTitle("CC #nu_{#mu} muon momentum fractional error (comparing recosim & truth)");
histmomsyst3->SetXTitle("P_{#mu,true}");
histmomsyst3->SetYTitle("#DeltaP_{#mu}/P_{#mu,true}");
histmomsyst3->GetYaxis()->SetTitleOffset(1.27);
histmomsyst3->GetYaxis()->SetLabelSize(.035);
set_plot_style();
// histmomsyst3->SetMinimum(histmomsyst3->GetMaximum()*.0);
histmomsyst3->GetXaxis()->SetTitleOffset(1.9);
histmomsyst3->GetYaxis()->SetTitleOffset(1.9);
histmomsyst3->Draw("lego2");
gPad->GetCanvas()->Print("histmomsyst2.pdf");
  gStyle->SetOptStat(0);
new TCanvas;

histmomsyst3->SetTitle("CC #nu_{#mu} muon momentum fractional error (comparing recosim & truth)");
histmomsyst3->SetXTitle("P_{#mu,true}");
histmomsyst3->SetYTitle("#DeltaP/P_{#mu,true}");
//histmomsyst3->GetYaxis()->SetTitleOffset(1.27);
histmomsyst3->GetYaxis()->SetLabelSize(.035);
set_plot_style();
histmomsyst3->GetXaxis()->SetRangeUser(5,25);
histmomsyst3->Draw("lego2");
gPad->GetCanvas()->Print("histmomsyst22.pdf");


new TCanvas;
histthetasyst3->SetTitle("CC #nu_{#mu} muon #theta_{#mu} fractional error (comparing recosim & truth)");
histthetasyst3->SetXTitle("#theta_{true} (degrees)");
histthetasyst3->SetYTitle("#Delta#theta/#theta_{true}");
histthetasyst3->GetYaxis()->SetTitleOffset(1.27);
histthetasyst3->GetYaxis()->SetLabelSize(.035);
// histthetasyst3->SetMinimum(histthetasyst3->GetMaximum()*.0);
histthetasyst3->GetXaxis()->SetTitleOffset(1.9);
histthetasyst3->GetYaxis()->SetTitleOffset(1.9);
histthetasyst3->Draw("lego2");
gPad->GetCanvas()->Print("histthetasyst2.pdf");

new TCanvas;
histthetasyst3->SetTitle("CC #nu_{#mu} muon #theta_{#mu} fractional error (comparing recosim & truth)");
histthetasyst3->SetXTitle("#theta_{true} (degrees)");
histthetasyst3->SetYTitle("#Delta#theta/#theta_{true}");
histthetasyst3->GetYaxis()->SetLabelSize(.035);
// histthetasyst3->GetXaxis()->SetRangeUser(.8,.95);
histthetasyst3->Draw("lego2");
gPad->GetCanvas()->Print("histthetasyst22.pdf");



TH1D* histmomsystfinal = new TH1D("histmomsystfinal", "",40,0.,50);

TH1D *py = histmomsyst2ext->ProjectionY("pY", 1, 1);
for(int i=1;i<=40;i++)
{
 //new TCanvas;
py = histmomsyst2ext->ProjectionY("pY", i, i); // where firstYbin = 0 and lastYbin = 9

double minf=-1., maxf=1.;
TF1 *gfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
py->Fit("Gaussian","RQN"); 
Double_t sigma = gfit->GetParameter(2); //value of 1st parameter
Double_t esigma = gfit->GetParError(2); //error on 1st parameter
histmomsystfinal->Fill(histmomsystfinal->GetBinCenter(i),sigma);
histmomsystfinal->SetBinError(i,esigma);
// py->Draw();
}


TH1D* histthetasystfinal = new TH1D("histthetasystfinal", "",36,0.,72);

TH1D *pyy = histthetasyst2ext->ProjectionY("pYY", 1, 1);
for(int i=1;i<=histthetasyst2ext->GetNbinsX();i++)
{
// new TCanvas;
pyy = histthetasyst2ext->ProjectionY("pYY", i, i); // where firstYbin = 0 and lastYbin = 9

double minf=-1., maxf=1.;
TF1 *gfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
pyy->Fit("Gaussian","RQN"); 
Double_t sigma = gfit->GetParameter(2); //value of 1st parameter
Double_t esigma = gfit->GetParError(2); //error on 1st parameter
histthetasystfinal->Fill(histthetasystfinal->GetBinCenter(i),sigma);
std::cout<<"large "<<histthetasystfinal->GetBinCenter(i)<<" "<<sigma<<std::endl;
histthetasystfinal->SetBinError(i,esigma);
//pyy->Draw();
}



TH1D* histmomsystfinall = new TH1D("histmomsystfinall", "",20,0.,25);

TH1D *ppy = histmomsyst2->ProjectionY("ppY", 1, 1);
for(int i=1;i<=histmomsyst2->GetNbinsX();i++)
{
// new TCanvas;
ppy = histmomsyst2->ProjectionY("ppY", i, i); // where firstYbin = 0 and lastYbin = 9

double minf=-1., maxf=1.;
TF1 *ggfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
ppy->Fit("Gaussian","RQN"); 
Double_t sigma = ggfit->GetParameter(2); //value of 1st parameter
Double_t esigma = ggfit->GetParError(2); //error on 1st parameter
histmomsystfinall->Fill(histmomsystfinall->GetBinCenter(i),sigma);
histmomsystfinall->SetBinError(i,esigma);
 //ppy->Draw();
}


new TCanvas;
ppy = histmomsyst2->ProjectionY("ppY", 9, 9);
TF1 *ggfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
ppy->SetTitle("CC #nu_{#mu} P_{#mu} recosim and truth, 10-11.25 GeV/c");
ppy->SetXTitle("(#Delta P_{#mu})/P_{#mu,true}");
ppy->SetYTitle("arbitrary units");
ppy->GetXaxis()->SetRangeUser(-2,2);
ppy->GetYaxis()->SetTitleOffset(1.27);
ppy->Rebin(2);
ppy->Fit("Gaussian","RQ"); 
 gStyle->SetOptStat(1000001110);
gPad->GetCanvas()->Print("momres.pdf");

TH1D* histthetasystfinall = new TH1D("histthetasystfinall", "",18,0.,36);

TH1D *ppyy = histthetasyst2->ProjectionY("ppYY", 1, 1);
for(int i=1;i<=histthetasyst2->GetNbinsX();i++)
{
// new TCanvas;
ppyy = histthetasyst2->ProjectionY("ppYY", i, i); // where firstYbin = 0 and lastYbin = 9

double minf=-1., maxf=1.;
TF1 *ggfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
ppyy->SetTitle("testtitle");
ppyy->Fit("Gaussian","RQN"); 
Double_t sigma = ggfit->GetParameter(2); //value of 1st parameter
Double_t esigma = ggfit->GetParError(2); //error on 1st parameter
histthetasystfinall->Fill(histthetasystfinall->GetBinCenter(i),sigma);
histthetasystfinall->SetBinError(i,esigma);
std::cout<<"small "<<histthetasystfinall->GetBinCenter(i)<<" "<<sigma<<std::endl;
//ppyy->Draw();
}

new TCanvas;
ppyy = histthetasyst2->ProjectionY("ppYY", 15, 15);
TF1 *ggfit = new TF1("Gaussian","gaus",minf,maxf); // Create the fit function
ppyy->SetTitle("CC #nu_{#mu} #theta_{#mu} recosim and truth, 28-30^{0}");
ppyy->SetXTitle("(#Delta#theta_{#mu})/#theta_{#mu,true}");
ppyy->SetYTitle("arbitrary units");
ppyy->GetXaxis()->SetRangeUser(-2,2);
ppyy->GetYaxis()->SetTitleOffset(1.27);
ppyy->Rebin(1);
ppyy->Fit("Gaussian","RQ"); 
gPad->GetCanvas()->Print("thetares.pdf");

new TCanvas;
 gStyle->SetOptStat(0);
histmomsystfinall->SetTitle("CC #nu_{#mu} muon momentum fractional error (comparing recosim & truth)");
histmomsystfinall->SetXTitle("P_{#mu,true}");
histmomsystfinall->SetYTitle("#delta[|#DeltaP_{#mu}|/P_{#mu,true}]");
histmomsystfinall->GetYaxis()->SetTitleOffset(1.27);
histmomsystfinall->GetYaxis()->SetLabelSize(.035);
histmomsystfinall->SetMinimum(0);
histmomsystfinall->Draw("E1");
gPad->GetCanvas()->Print("profmomsyst.pdf");
new TCanvas;

   histthetasystfinall->SetTitle("CC #nu_{#mu} muon #theta_{#mu} fractional error (comparing recosim & truth)");
   histthetasystfinall->SetXTitle("#theta_{#mu, true}");
   histthetasystfinall->SetYTitle("#delta[|#Delta#theta_{#mu}|/#theta_{#mu, true}]");
   histthetasystfinall->GetYaxis()->SetTitleOffset(1.27);
   histthetasystfinall->GetYaxis()->SetLabelSize(.035);
   histthetasystfinall->SetMinimum(0);
histthetasystfinall->Draw("E1");
gPad->GetCanvas()->Print("profthetasyst.pdf");

TFile *g=new TFile("profsyst.root", "RECREATE");
histmomsystfinal->Write();
histthetasystfinal->Write();
g->Close();


// new TCanvas; 
// histthetasystfinal->GetXaxis()->SetRangeUser(.85,1.);
// histthetasystfinal->GetYaxis()->SetLabelSize(.035);
// histthetasystfinal->DrawCopy("E1");
// gPad->GetCanvas()->Print("profthetasyst2.pdf");


//   new TCanvas;
//   histtest->DrawCopy();
//  
//   new TCanvas;
//   histtest2->DrawCopy(); 
//   
//   new TCanvas;
//   histtest3->DrawCopy(); 
//   
//     new TCanvas;
//   histtest4->DrawCopy("colz");

new TCanvas;
study->SetTitle("CC");
study->Draw("colz");
std::cout<<study->GetMaximum()*.005<<std::endl;



new TCanvas;
study_qe->SetTitle("CCQE");
study_qe->Draw("colz");
std::cout<<study->GetMaximum()*.005<<std::endl;


new TCanvas;
study_res->SetTitle("CCRES");
study_res->Draw("colz");
std::cout<<study->GetMaximum()*.005<<std::endl;


new TCanvas;
study_dis->SetTitle("CCDIS");
study_dis->Draw("colz");
std::cout<<study->GetMaximum()*.005<<std::endl;


new TCanvas;

hist4->Divide(hist1);
hist5->Divide(hist1);
hist6->Divide(hist1);
hist4->SetMinimum(0);
 hist4->DrawNormalized("", hist4->Integral()*potfactor);
 hist5->SetLineColor(7);
 hist5->DrawNormalized("SAME", hist5->Integral()*potfactor);
 hist6->SetLineColor(6);
 hist6->DrawNormalized("SAME", hist6->Integral()*potfactor);


new TCanvas;

histeff_qe2->Divide(histeff_qe1);
histeff_res2->Divide(histeff_res1);
histeff_dis2->Divide(histeff_dis1);

histeff_qe2->SetLineColor(3);
histeff_qe2->SetMinimum(0);
histeff_qe2->SetMaximum(1);
histeff_qe2->Draw();
histeff_res2->SetLineColor(7);
histeff_res2->Draw("SAME");
histeff_dis2->SetLineColor(6);
histeff_dis2->Draw("SAME");

new TCanvas;
hit2->SetLineColor(1);
hit2->DrawNormalized("", hit2->Integral()*potfactor);
hit2iso->SetLineColor(2);
hit2iso->DrawNormalized("SAME", hit2iso->Integral()*potfactor);
new TCanvas;
hit3->SetLineColor(1);
hit3->DrawNormalized("", hit3->Integral()*potfactor);
hit3iso->SetLineColor(2);
hit3iso->DrawNormalized("SAME", hit3iso->Integral()*potfactor);




new TCanvas;

hit4->SetLineColor(1);
hit4->DrawNormalized("", hit4->Integral()*potfactor);
hit4iso->SetLineColor(2);
hit4iso->DrawNormalized("SAME", hit4iso->Integral()*potfactor);
// TFile *momcosiso=new TFile("momcosiso.root", "RECREATE");
new TCanvas;

hit3iso->SetTitle("Argon to isoscalar correction factor");
hit3iso->SetXTitle("P_{#mu} (GeV/c)");
hit3iso->SetYTitle("Isoscalar correction factor");
hit3iso->Divide(hit3);
hit3iso->SetLineColor(1);
hit3iso->GetYaxis()->SetTitleOffset(1.3);
hit3iso->SetMaximum(1);
hit3iso->SetMinimum(0.9);
 TF1 *fitit = new TF1("","pol1",0,25); 
fitit->SetParNames("intercept","slope");
hit3iso->Fit(fitit,"EMR");
hit3iso->Draw("hist");
fitit->Draw("SAME");
// hit3iso->Write();
for(int i=1;i<=hit3iso->GetNbinsX();i++)
std::cout<<"mom "<<hit3iso->GetBinCenter(i)<<" "<<hit3iso->GetBinContent(i)<<std::endl;
//  gPad->GetCanvas()->Print("isomom.pdf");
new TCanvas;
hit4iso->SetTitle("Argon to isoscalar correction factor");
hit4iso->SetXTitle("#theta_{#mu} (degrees)");
hit4iso->SetYTitle("Isoscalar correction factor");
hit4iso->Divide(hit4);
hit4iso->SetLineColor(1);
hit4iso->GetYaxis()->SetTitleOffset(1.3);
hit4iso->SetMaximum(1);
hit4iso->SetMinimum(0.9);
 TF1 *fitit2 = new TF1("","pol1",0,40); 
fitit2->SetParNames("intercept","slope");
hit4iso->Fit(fitit2,"EMR");
hit4iso->Draw("hist");
fitit2->Draw("SAME");
// hit4iso->Write();
// momcosiso->Close();
//  gPad->GetCanvas()->Print("isocos.pdf");
for(int i=1;i<=hit4iso->GetNbinsX();i++)
std::cout<<"cos "<<hit4iso->GetBinCenter(i)<<" "<<hit4iso->GetBinContent(i)<<std::endl;
new TCanvas;
hit2iso->SetTitle("Argon to isoscalar correction factor");
hit2iso->SetXTitle("E_{#nu} (GeV)");
hit2iso->SetYTitle("Isoscalar correction factor");
hit2iso->Divide(hit2);
hit2iso->SetLineColor(1);
hit2iso->GetYaxis()->SetTitleOffset(1.3);
hit2iso->SetMaximum(1);
hit2iso->SetMinimum(0.9);
 TF1 *fitit3 = new TF1("","pol1",0,50); 
fitit3->SetParNames("intercept","slope");
hit2iso->Fit(fitit3,"EMR");
hit2iso->Draw("hist");
fitit3->Draw("SAME");
//  gPad->GetCanvas()->Print("isoenergy.pdf");
for(int i=1;i<=hit2iso->GetNbinsX();i++)
std::cout<<"energy "<<hit2iso->GetBinCenter(i)<<" "<<hit2iso->GetBinContent(i)<<std::endl;

 new TCanvas;

 histunfoldmom2->SetTitle("CC #nu_{#mu} muon momentum recosim - truth (after cuts)");
 gStyle->SetOptStat("mr");
gPad->GetCanvas()->Update();
 histunfoldmom2->SetXTitle("P_{#mu,recosim} - P_{#mu,true} (GeV/c)");
 histunfoldmom2->GetYaxis()->SetTitleOffset(1.27);
 histunfoldmom2->SetYTitle("fraction of events");
 histunfoldmom2->Scale(1./histunfoldmom2->Integral());
 histunfoldmom22->Scale(1./histunfoldmom22->Integral());
//  histunfoldmom2->Fit("gaus","W","",-2.5,2.5);

histunfoldmom2->Draw();
gPad->Modified();
gPad->Update();
TPaveStats *st1 = (TPaveStats*)histunfoldmom2->GetListOfFunctions()->FindObject("stats");
st1->SetLineColor(1);
st1->SetY1NDC(0.84);
st1->SetY2NDC(0.98);
histunfoldmom22->SetLineColor(2);


histunfoldmom22->Draw("SAMES");
gPad->Modified();
gPad->Update();
TPaveStats *st2 = (TPaveStats*)histunfoldmom22->GetListOfFunctions()->FindObject("stats");
st2->SetLineColor(2);
st2->SetY1NDC(0.68);
st2->SetY2NDC(0.82);
gPad->Modified();

laex->Draw();
laexx->Draw();

 gPad->GetCanvas()->Print("histunfoldmom2.pdf");
 


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


void
set_plot_style()
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
